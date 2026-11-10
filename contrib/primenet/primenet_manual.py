#!/usr/bin/env python3
"""Automate the PrimeNet *manual* testing workflow (assignment + result upload).

This drives the same web pages a logged-in user would use, so it does not
impersonate a trusted client: source-built mprime cannot use the automatic
PrimeNet protocol (server replies "pnErrorResult=9 Untrusted program versions
currently excluded by PrimeNet"), but the manual-testing forms are open to any
account.

Usage:
    set PN_PASS=...                 # PrimeNet password (never written to disk)
    python primenet_manual.py login --user vm001
    python primenet_manual.py assign --pref 5 --exp-lo 1277 --exp-hi 1277 \
        --out worktodo.txt
    python primenet_manual.py submit --file results.json.txt

Cookie jar is stored next to this script as .pn_cookies.txt.
"""
import argparse, os, re, sys, urllib.parse, urllib.request, http.cookiejar

BASE = "https://www.mersenne.org"
UA = "Mozilla/5.0 (mprime manual-testing helper)"
STATE_DIR = os.environ.get("PN_STATE", os.path.dirname(os.path.abspath(__file__)))
COOKIES = os.path.join(STATE_DIR, ".pn_cookies.txt")
PREF_NAMES = {
    150: "first-time PRP", 151: "PRP double-check", 152: "world-record PRP",
    153: "100M-digit PRP", 155: "PRP double-check with proof",
    4: "P-1 factoring", 5: "ECM factoring", 2: "trial factoring",
    8: "ECM on cofactors", 160: "cofactor first PRP", 161: "cofactor PRP DC",
    101: "LL double-check",
}
ASSIGN_RE = re.compile(r"(?:^|>)\s*((?:Test|DoubleCheck|PRPDC|PRP|Pminus1|Pfactor|ECM2?|Factor|P\+1|CERT)=[^\r\n<]+)",
                       re.M)


def make_opener():
    jar = http.cookiejar.MozillaCookieJar(COOKIES)
    if os.path.exists(COOKIES):
        try:
            jar.load(ignore_discard=True, ignore_expires=True)
        except Exception:
            pass
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
    opener.addheaders = [("User-Agent", UA)]
    return opener, jar


def save_jar(jar):
    try:
        jar.save(ignore_discard=True, ignore_expires=True)
    except Exception:
        pass


def post(opener, path, fields):
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(BASE + path, data=data)
    req.add_header("Content-Type", "application/x-www-form-urlencoded")
    return opener.open(req, timeout=60).read().decode("utf-8", "replace")


def logged_in_html(html):
    return "logged in" in html or "Logout" in html


def cmd_login(opener, jar, user, password):
    html = opener.open(BASE + "/update/", timeout=60).read().decode("utf-8", "replace")
    html = post(opener, "/update/", {"user_login": user, "user_password": password})
    if logged_in_html(html):
        m = re.search(r'name="uid" value="([0-9]+)"', html)
        save_jar(jar)
        print("login ok" + (" (uid %s)" % m.group(1) if m else ""))
        return True
    print("login FAILED")
    return False


def cmd_assign(opener, jar, pref, workers, num, lo, hi, out, append=True):
    fields = {"workers": str(workers), "num_to_get": str(num), "pref": str(pref),
              "exp_lo": str(lo or ""), "exp_hi": str(hi or "")}
    html = post(opener, "/manual_assignment/", fields)
    lines = ASSIGN_RE.findall(html)
    if not lines:
        txt = re.sub(r"<[^>]+>", " ", html)
        txt = re.sub(r"\s+", " ", txt)
        i = txt.find("PrimeNet Get Manual Assignments")
        print("no assignments returned:", txt[i:i + 300] if i >= 0 else txt[:300])
        return 2
    print("\n".join(lines))
    if out:
        mode = "a" if append and os.path.exists(out) else "w"
        with open(out, mode, encoding="utf-8") as f:
            for ln in lines:
                f.write(ln + "\n")
        print("appended %d assignment(s) to %s" % (len(lines), out))
    return 0


def multipart(fields):
    bnd = "----primenetmanual7f3a91"
    out = []
    for k, v in fields.items():
        out.append("--" + bnd)
        if k == "data_file":
            out.append('Content-Disposition: form-data; name="data_file"; filename=""')
            out.append("Content-Type: application/octet-stream")
            out.append("")
            out.append("")
        else:
            out.append('Content-Disposition: form-data; name="%s"' % k)
            out.append("")
            out.append(v)
    out.append("--" + bnd + "--")
    out.append("")
    return ("\r\n".join(out)).encode("utf-8"), bnd


def cmd_submit(opener, jar, path):
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    html = opener.open(BASE + "/manual_result/", timeout=60).read().decode("utf-8", "replace")
    m = re.search(r'name="was_logged_in_as" value="([0-9]+)"', html)
    if not m:
        print("not logged in / no login field on manual_result page")
        return 2
    fields = {"was_logged_in_as": m.group(1), "data": text, "data_file": ""}
    body, bnd = multipart(fields)
    req = urllib.request.Request(BASE + "/manual_result/", data=body)
    req.add_header("Content-Type", "multipart/form-data; boundary=" + bnd)
    resp = opener.open(req, timeout=120).read().decode("utf-8", "replace")
    try:
        with open(os.path.join(STATE_DIR, "last_submit.html"), "w", encoding="utf-8") as f:
            f.write(resp)
    except OSError:
        pass
    plain = re.sub(r"<[^>]+>", " ", resp)
    plain = re.sub(r"\s+", " ", plain)
    keys = ("accepted", "Accepted", "rejected", "Rejected", "Error", "error", "GHz", "processed", "verified", "Credited")
    hits = [plain.find(k) for k in keys if plain.find(k) >= 0]
    i = min(hits) if hits else 0
    print(plain[i:i + 800])
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--user", default=os.environ.get("PN_USER", "vm001"))
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("login"); p.add_argument("--password", default=os.environ.get("PN_PASS"))
    p = sub.add_parser("assign")
    p.add_argument("--pref", type=int, default=155)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--num", type=int, default=1)
    p.add_argument("--exp-lo", type=int); p.add_argument("--exp-hi", type=int)
    p.add_argument("--out", default="worktodo.txt")
    p = sub.add_parser("submit"); p.add_argument("--file", required=True)
    a = ap.parse_args()
    opener, jar = make_opener()
    if a.cmd == "login":
        if not a.password:
            print("password required (--password or PN_PASS)"); return 2
        return 0 if cmd_login(opener, jar, a.user, a.password) else 1
    if a.cmd == "assign":
        return cmd_assign(opener, jar, a.pref, a.workers, a.num, a.exp_lo, a.exp_hi, a.out)
    if a.cmd == "submit":
        return cmd_submit(opener, jar, a.file)
    return 2


if __name__ == "__main__":
    sys.exit(main())
