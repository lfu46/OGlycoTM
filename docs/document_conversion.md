# Document conversion recipes

> Moved out of `CLAUDE.md` on 2026-09-02 so it is loaded on demand rather than into every session.

The Markdown-to-PDF recipe (xelatex, Unicode caveats, write to local disk first) also lives in the global `~/.claude/CLAUDE.md`; this file keeps the project-specific variants (Times New Roman + TOC, the subscript sed filter) and the DOCX, TIFF and text recipes.

## Markdown to PDF Conversion

Convert markdown files to text-selectable PDF using pandoc:

```bash
pandoc "input.md" -o "output.pdf" --pdf-engine=xelatex \
  -V geometry:margin=1in -V fontsize=11pt -V mainfont="Times New Roman" \
  --toc -V colorlinks=true -V linkcolor=blue -V urlcolor=blue
```

If Unicode subscripts/superscripts cause issues (Times New Roman doesn't support them), pipe through sed first:
```bash
sed -e 's/⁻/-/g' -e 's/₀/0/g' -e 's/₁/1/g' -e 's/₂/2/g' -e 's/₃/3/g' -e 's/₄/4/g' -e 's/₅/5/g' -e 's/₆/6/g' -e 's/₇/7/g' -e 's/₈/8/g' -e 's/₉/9/g' "input.md" | \
  pandoc -o "output.pdf" --pdf-engine=xelatex \
  -V geometry:margin=1in -V fontsize=11pt -V mainfont="Times New Roman" \
  --toc -V colorlinks=true -V linkcolor=blue -V urlcolor=blue
```

## DOCX to PDF Conversion

Convert Word documents to text-selectable PDF using pandoc:

```bash
pandoc "input.docx" -o "output.pdf" --pdf-engine=xelatex
```

Note: This produces text-selectable PDFs but may not preserve complex formatting (tables, images, custom styles) perfectly. For exact formatting preservation, use `docx2pdf` (`pip install docx2pdf`) which requires Microsoft Word or LibreOffice.

## PDF to High-Resolution TIFF Conversion

Convert PDF files to high-resolution TIFF (600 DPI) using ImageMagick:

```bash
# Single file
magick -density 600 "input.pdf" -quality 100 "output.tiff"

# Batch convert all PDFs in a folder
for pdf in /path/to/folder/*.pdf; do
  filename=$(basename "$pdf" .pdf)
  magick -density 600 "$pdf" -quality 100 "/path/to/output/${filename}.tiff"
done
```

Note: Requires ImageMagick (`brew install imagemagick`). 600 DPI TIFFs are large (~80-90 MB each) but publication-quality.

## PDF to Text Extraction (Token-Efficient)

Extract text from PDF files for uploading to Claude Desktop (avoids image tokens):

```bash
# Requires poppler: brew install poppler

# Single file (compact, token-efficient)
pdftotext "input.pdf" "output.txt"

# Batch convert all PDFs in a folder
for pdf in /path/to/folder/*.pdf; do
  pdftotext "$pdf" "${pdf%.pdf}.txt"
done
```

Note: The default mode (no flags) is most token-efficient. Use `-layout` only if you need to preserve table/column structure (costs more tokens due to extra whitespace).
