

### Reading the htseq count files into R was having an issue related to columns in the counts file and not matching, to fix this, just duplicate the first column (gene name). This code should work recursively within a directory over each tissue type. Substitute lepSecond_counts_files with bryFirst_counts_files to modify those

```
find lepSecond_counts_files -type f -name '*.tsv' -exec zsh -c '
    for file do
        awk -F"\t" '\''BEGIN {OFS = FS} $2 == "" { $2 = $1 } { print $0 }'\'' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
    done
' zsh {} +
```
