#!/usr/bin/env python3
"""Top up the farm worktodo when any [Worker #n] section runs low.
Fetches ECM assignments via the PrimeNet manual form, expands each into
single-curve ECM2 lines (results are batched per line by mprime), rewrites
worktodo.txt and lets the watchdog restart mprime (workers re-read the file
at startup; killed workers resume from their save files).
Runs ON the box. Stdlib only."""
import os, re, sys, urllib.parse, urllib.request, http.cookiejar

BASE = "https://www.mersenne.org"
HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.dirname(HERE)
FARM = os.path.join(MP, "farm1")
WT = os.path.join(FARM, "worktodo.txt")
COOKIES = os.path.join(HERE, ".pn_cookies.txt")
PASSFILE = os.path.join(MP, ".pn_pass")
USER = os.environ.get("PN_USER", "vm001")
MIN_LINES = 3          # top up when a section drops below this
FETCH = 20             # assignments to fetch per top-up
EXP_LO, EXP_HI = 100000, 400000
ASSIGN_RE = re.compile(r"(?:^|>)\s*((?:Test|DoubleCheck|PRPDC|PRP|Pminus1|Pfactor|ECM2?|Factor|P\+1|CERT)=[^\r\n<]+)", re.M)

def password():
    if os.environ.get("PN_PASS"):
        return os.environ["PN_PASS"]
    if os.path.exists(PASSFILE):
        return open(PASSFILE).read().strip()
    return None

def section_counts(lines):
    counts, cur = {}, None
    for l in lines:
        m = re.match(r"\[(Worker #\d+)\]", l)
        if m:
            cur = m.group(1)
            counts[cur] = 0
        elif cur and l.strip().startswith(("ECM2=", "ECM=", "Pminus1=", "Pfactor=", "PRP", "Test=", "DoubleCheck=", "Factor=", "CERT=")):
            counts[cur] += 1
    return counts

def split_assignment(line):
    """ECM2=k,b,n,c,B1,B2,curves[,...] -> list of single-curve lines."""
    if not line.startswith("ECM2="):
        return [line]
    parts = line.split(",")
    if len(parts) < 8:
        return [line]
    try:
        curves = int(parts[7])
    except ValueError:
        return [line]
    if curves <= 1:
        return [line]
    return [",".join(parts[:7] + ["1"]) for _ in range(curves)]

def main():
    lines = open(WT, encoding="utf-8", errors="replace").read().splitlines()
    counts = section_counts(lines)
    low = [s for s, c in counts.items() if c < MIN_LINES]
    if not low:
        print("sections healthy:", counts)
        return 0
    print("low sections:", low, counts)

    jar = http.cookiejar.MozillaCookieJar(COOKIES)
    if os.path.exists(COOKIES):
        try:
            jar.load(ignore_discard=True, ignore_expires=True)
        except Exception:
            pass
    op = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
    op.addheaders = [("User-Agent", "mprime-arm64-farm/1.0")]

    def post(path, fields):
        req = urllib.request.Request(BASE + path,
                                     data=urllib.parse.urlencode(fields).encode())
        req.add_header("Content-Type", "application/x-www-form-urlencoded")
        return op.open(req, timeout=120).read().decode("utf-8", "replace")

    html = post("/manual_assignment/", {"workers": "1", "num_to_get": str(FETCH),
                                        "pref": "5", "exp_lo": str(EXP_LO),
                                        "exp_hi": str(EXP_HI)})
    got = ASSIGN_RE.findall(html)
    if not got:
        m = re.search(r"name=\"was_logged_in_as\" value=", html)
        if not m:
            pw = password()
            if not pw:
                print("not logged in and no password")
                return 2
            post("/update/", {"user_login": USER, "user_password": pw})
            try:
                jar.save(ignore_discard=True, ignore_expires=True)
            except Exception:
                pass
            html = post("/manual_assignment/", {"workers": "1", "num_to_get": str(FETCH),
                                                "pref": "5", "exp_lo": str(EXP_LO),
                                                "exp_hi": str(EXP_HI)})
            got = ASSIGN_RE.findall(html)
    if not got:
        print("no assignments returned")
        return 1

    newlines = []
    for a in got:
        newlines.extend(split_assignment(a.strip()))
    print("fetched %d assignments -> %d lines" % (len(got), len(newlines)))

    # distribute: append round-robin to the low sections
    li = 0
    out = []
    for l in lines:
        out.append(l)
        m = re.match(r"\[(Worker #\d+)\]", l)
        if m:
            sec = m.group(1)
            if sec in low:
                # insert at end of this section = before next header; handled by marker
                pass
    # simpler rebuild: group sections
    sections, order, cur = {}, [], None
    for l in lines:
        m = re.match(r"\[(Worker #\d+)\]", l)
        if m:
            cur = m.group(1)
            order.append(cur)
            sections[cur] = []
        elif cur:
            sections[cur].append(l)
        else:
            order.append(("__root__", l)) if False else None
    i = 0
    for sec in order:
        take = [nl for j, nl in enumerate(newlines) if j % max(1, len(low)) == i]
        sections[sec].extend(take)
        i += 1
    final = []
    for sec in order:
        final.append("[%s]" % sec)
        final.extend([l for l in sections[sec] if l.strip()])
    tmp = WT + ".new"
    with open(tmp, "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(final) + "\n")
    os.replace(tmp, WT)
    print("worktodo updated; killing mprime so it re-reads (watchdog will restart)")
    os.system("pkill -x mprime")
    return 0

if __name__ == "__main__":
    sys.exit(main())
