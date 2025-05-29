dir="inputs/scRNA-seq"
for type in `ls $dir`; do
    for celltype in `ls $dir/$type`; do
        for folder in `ls $dir/$type/$celltype`; do
            if [[ -d $dir/$type/$celltype/$folder ]]; then
                if [[ "$folder" =~ [\D\W]*TFs-500[\D\W]* ]]; then
                    echo "This is a 500 folder: $dir/$type/$celltype/$folder"
                    mv $dir/$type/$celltype/$folder $dir/$type/$celltype/TFs-500
                else
                    echo "This is a 1000 folder: $dir/$type/$celltype/$folder"
                    mv $dir/$type/$celltype/$folder $dir/$type/$celltype/TFs-1000
                fi
            else
                echo ""
            fi
        done
    done
done