#!/bin/bash
for f in *.csv; do
    tmp="${f}.tmp"
    (echo -ne '\xEF\xBB\xBF' && iconv -f euc-kr -t utf-8 -c "$f") > "$tmp" && mv "$tmp" "$f"
done

