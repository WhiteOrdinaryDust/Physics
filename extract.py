import sys
from pypdf import PdfReader
r = PdfReader(r"D:\DeepSeek Harness\observer-scan\thooft_cai.pdf")
print("pages:", len(r.pages))
out = []
for i, p in enumerate(r.pages):
    t = p.extract_text() or ""
    out.append("\n=== PAGE %d ===\n" % (i+1) + t)
txt = "".join(out)
open(r"D:\DeepSeek Harness\observer-scan\thooft_cai.txt","w",encoding="utf-8").write(txt)
print("chars:", len(txt))
