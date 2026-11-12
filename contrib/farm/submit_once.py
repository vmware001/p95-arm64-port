#!/usr/bin/env python3
"""Submit farm results.json.txt to PrimeNet (manual workflow), adding the
checksum object the /manual_result/ form requires since 2026.
Runs ON the box. Password: /data/adb/mprime/.pn_pass (root 600) or PN_PASS.
Cookies: <script_dir>/.pn_cookies.txt. Stdlib only."""
import json, os, re, sys, urllib.parse, urllib.request, http.cookiejar, zlib

BASE = "https://www.mersenne.org"
HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.dirname(HERE)
FARM = os.path.join(MP, "farm1")
COOKIES = os.path.join(HERE, ".pn_cookies.txt")
PASSFILE = os.path.join(MP, ".pn_pass")
RESULTS = os.path.join(FARM, "results.json.txt")
USER = os.environ.get("PN_USER", "vm001")
UA = "mprime-arm64-farm/1.0"

LAYOUT = "{exp};ECM;;{curves};{b1};{b2};;;;{prog};{ver};;;ARM64 Linux;{ts}"

def password():
    if os.environ.get("PN_PASS"):
        return os.environ["PN_PASS"]
    if os.path.exists(PASSFILE):
        return open(PASSFILE).read().strip()
    return None

def add_checksum(line):
    d = json.loads(line)
    s = LAYOUT.format(exp=d["exponent"], curves=d.get("curves", 1), b1=d["b1"],
                      b2=d["b2"], prog=d["program"]["name"],
                      ver=d["program"]["version"], ts=d["timestamp"])
    d["checksum"] = {"version": 1,
                     "checksum": "%08X" % (zlib.crc32(s.encode()) & 0xFFFFFFFF)}
    return json.dumps(d, separators=(", ", ": "))

def main():
    lines = [l.strip() for l in open(RESULTS, encoding="utf-8", errors="replace") if l.strip()]
    if not lines:
        print("no results")
        return 0
    text = "\n".join(add_checksum(l) for l in lines) + "\n"

    jar = http.cookiejar.MozillaCookieJar(COOKIES)
    if os.path.exists(COOKIES):
        try:
            jar.load(ignore_discard=True, ignore_expires=True)
        except Exception:
            pass
    op = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
    op.addheaders = [("User-Agent", UA)]

    def get(path):
        return op.open(BASE + path, timeout=90).read().decode("utf-8", "replace")

    def post(path, fields, multipart=False):
        if multipart:
            bnd = "----p95farm"
            body = b""
            for k, v in fields.items():
                body += ("--%s\r\nContent-Disposition: form-data; name=\"%s\"\r\n\r\n%s\r\n"
                         % (bnd, k, v)).encode()
            body += ("--%s--\r\n" % bnd).encode()
            req = urllib.request.Request(BASE + path, data=body)
            req.add_header("Content-Type", "multipart/form-data; boundary=" + bnd)
        else:
            req = urllib.request.Request(BASE + path,
                                         data=urllib.parse.urlencode(fields).encode())
            req.add_header("Content-Type", "application/x-www-form-urlencoded")
        return op.open(req, timeout=120).read().decode("utf-8", "replace")

    html = get("/manual_result/")
    m = re.search(r'name="was_logged_in_as" value="([0-9]+)"', html)
    if not m:
        pw = password()
        if not pw:
            print("not logged in and no password available")
            return 2
        get("/update/")  # session warmup, mirrors _pnm.py
        r = post("/update/", {"user_login": USER, "user_password": pw})
        try:
            jar.save(ignore_discard=True, ignore_expires=True)
        except Exception:
            pass
        print("login page has logout:", "Logout" in r or "logged in" in r)
        print("cookies:", [(c.name, c.domain) for c in jar])
        try:
            jar.save(ignore_discard=True, ignore_expires=True)
        except Exception:
            pass
        # reload the jar cross-process style (proven path in _pnm.py)
        jar = http.cookiejar.MozillaCookieJar(COOKIES)
        jar.load(ignore_discard=True, ignore_expires=True)
        op = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
        op.addheaders = [("User-Agent", UA)]
        html = get("/manual_result/")
        m = re.search(r'name="was_logged_in_as" value="([0-9]+)"', html)
        if not m:
            print("login failed; server said:", re.sub(r"\s+", " ", r)[:200])
            print("manual_result now shows:", re.sub(r"\s+", " ", html)[:400])
            return 2

    resp = post("/manual_result/",
                {"was_logged_in_as": m.group(1), "data": text, "data_file": ""},
                multipart=True)
    plain = re.sub(r"<[^>]+>", " ", resp)
    plain = re.sub(r"\s+", " ", plain)
    i = plain.find("processing:")
    print(plain[i:i + 900] if i >= 0 else plain[:600])
    ok = ("Accepted" in plain and re.search(r"Accepted [1-9]", plain)) or \
         ("Done processing" in plain and "did not understand" not in plain.lower())
    print("SUBMIT-OK" if ok else "SUBMIT-UNCERTAIN")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
