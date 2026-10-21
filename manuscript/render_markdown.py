"""Create a readable Markdown companion from the authoritative LaTeX draft.

Requires Pandoc on PATH. Run after editing main.tex or references.bib.
"""
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
result = subprocess.run(
    ["pandoc", "main.tex", "--from=latex", "--to=gfm", "--standalone",
     "--template=markdown.template", "--citeproc", "--bibliography=references.bib",
     "--wrap=none"], cwd=HERE, capture_output=True, check=True, encoding="utf-8")
text = result.stdout
text = re.sub(r'<a href="#[^"]*"[^>]*>(.*?)</a>', r'\1', text)
text = re.sub(r'<figure[^>]*>\s*<embed src="([^"]+)"\s*/>\s*<figcaption>(.*?)</figcaption>\s*</figure>',
              lambda m: f"![{m.group(2).strip()}]({m.group(1).replace('.pdf', '.png')})",
              text, flags=re.S)
text = re.sub(r'<div id="refs"[^>]*>', '\n## References\n', text)
text = re.sub(r'</?div[^>]*>', '', text)
text = re.sub(r'\s*\{#tab:[^}]*\}', '', text)
text = re.sub(r'\\label\{[^}]*\}', '', text)
text = re.sub(r'``` math\n(.*?)\n```', r'$$\n\1\n$$', text, flags=re.S)
text = re.sub(r'\$`(.*?)`\$', r'$\1$', text)
text = re.sub(r'\n{3,}', '\n\n', text)
(HERE / 'article.md').write_text(text, encoding='utf-8')
print('Updated article.md from main.tex and references.bib.')
