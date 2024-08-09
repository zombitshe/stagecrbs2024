#!/usr/bin/env bash

input_file="$1" #input .b2g file path
output_file="$2" #give output file a name

touch "$output_file"
tmp_file=$(mktemp)

# Read every line
while IFS= read -r line || [[ -n "$line" ]]; do
    # Remove carriage return character at the end of $line (due to b2g line)
    line="${line%$'\r'}"
    
    # Check if the line contains RepID
    if grep -q 'RepID' <<< "$line"; then
        id=$(grep -oP '(?<=RepID=)(A0A[^_\s]*|[^\s]+)' <<< "$line") # For uncurated uniprot entries (start w/ A0A), need to stop before "_", however not if a curated uniprot entry 
        
        # Download JSON based on RepID type
        if [[ $id == UPI* ]]; then 
            url="https://rest.uniprot.org/uniparc/${id}.json"
            type="RefSeq"
        else
            url="https://rest.uniprot.org/uniprotkb/${id}.json"
            type="UniProt"
        fi
        
        # Download JSON and process data
        wget -q "$url" -O "${id}.json"
        echo "Proceeding .json for RepID=${id}"
        
        # Check if the file was downloaded successfully
        if [[ -f "${id}.json" ]]; then
            # if Refseq ID
            if [[ $type == "RefSeq" ]]; then 
            gene_name=$(jq -r '.uniParcCrossReferences[0].geneName' "${id}.json")
            access=$(jq -r '.uniParcCrossReferences[0].id' "${id}.json")
    
            printf "%s\tgene_name=%s\t${type}_access=%s\n" "$line" "$gene_name" "$access" >> "$tmp_file"
            rm "${id}.json"    # Code to execute if type is "RefSeq"
            else
            # if UniProt ID
            gene_name=$(jq -r '.genes[0].geneName.value' "${id}.json")
            access=$(jq -r '.primaryAccession' "${id}.json")
            
            printf "%s\tgene_name=%s\t${type}_access=%s\n" "$line" "$gene_name" "$access" >> "$tmp_file"
            rm "${id}.json"    # Code to execute if type is "RefSeq"
            fi
        else
            # If JSON download fails, just append the original line
            echo "$line" >> "$tmp_file"
        fi
    else
        # If line does not contain RepID, just append the original line
        echo "$line" >> "$tmp_file"
    fi
done < "$input_file"

# Move the temporary file to the output file
mv "$tmp_file" "$output_file"