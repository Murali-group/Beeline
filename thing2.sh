dir="config-files/scRNA-seq"
for celltype in `ls $dir`; do
    for config in `ls $dir/$celltype`; do
        # echo $dir/$celltype/$config
        if [[ -f $dir/$celltype/$config ]]; then
            echo "Now running configuration file: $dir/$celltype/$config";
            python BLEvaluator.py --config $dir/$celltype/$config -e & 
        else
            echo "Now skipping directory: $dir/$celltype/$config";
        fi
    done
done