

##reading the htseq count files into R was having an issue related to columns in the counts file and not matching, to fix this, just duplicate the first column (gene name), and it should work

```
for file in *.tsv; do    
    awk -F'\t' 'BEGIN {OFS = FS} $2 == "" { $2 = $1 } { print $0 }' "$file" > "$file.tmp"
    mv "$file.tmp" "$file"
done
```
