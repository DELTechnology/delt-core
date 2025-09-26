pandoc protocols.md -o paper.pdf \
  --pdf-engine=xelatex \
  -V geometry:"top=1.2in, bottom=1.2in, left=1.25in, right=1.25in" \
  -V fontsize=11pt \
  -V linestretch=1.15 \
  --number-sections \
  -V colorlinks=true \
  -V linkcolor=blue \
  -V urlcolor=blue

pandoc protocols.md -o protocols.docx