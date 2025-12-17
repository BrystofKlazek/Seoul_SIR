#!/bin/bash
for f in data_*시.csv; do
    new=$(echo "$f" | sed -E 's/시\.csv$/h.csv/')
    mv "$f" "$new"
done

