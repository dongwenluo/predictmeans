## RTF output window saved by Genstat -> plain text (the line with the working
## directory is dropped). Usage: python3 -I dev/genstat/rtf2txt.py <in.rtf> <out.txt>
import re, sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src, encoding="cp1252", errors="replace").read()
t = re.sub(r"\{\\fonttbl.*?\}\}", "", t, flags=re.S)
t = re.sub(r"\{\\colortbl.*?\}", "", t, flags=re.S)
t = re.sub(r"(\\[a-zA-Z]+-?\d*)\r?\n", r"\1 ", t)   # a line end delimits a control word
t = t.replace("\r", "").replace("\n", "")             # other RTF line breaks carry no meaning
t = t.replace("\\\\", "\x00")                         # literal backslash
for w, rep in (("par", "\n"), ("line", "\n"), ("page", "\n"), ("tab", "\t")):
    t = re.sub(r"\\" + w + r"(?![a-z])\s?", rep, t)
t = re.sub(r"\\[a-zA-Z]+-?\d*\s?", "", t)             # other control words
t = t.replace("{", "").replace("}", "").replace("\x00", "\\")
lines = [l.rstrip() for l in t.split("\n")]
lines = [l for l in lines if "WORKINGDIRECTORY" not in l]   # the user's folder
out, blank = [], 0
for l in lines:
    blank = blank + 1 if not l.strip() else 0
    if blank <= 1: out.append(l)
open(dst, "w").write("\n".join(out) + "\n")
