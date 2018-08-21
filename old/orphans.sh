
exit 0

# Identify orphans
for type in samples jobs ; do
    for dir in /mnt/wigclust*/data/unsafe/paa/mums-new/$type ; do
        for file in $(cd $dir ; rmdir */ 2> /dev/null ; ls) ; do            
            ourcopy=$(readlink $type/$file)
            if [ ! -e $type/$file ] || [ "$ourcopy" != $dir/$file ] ; then
                mkdir -p orphans
                ln -sfT $dir/$file orphans/$file
            fi
        done
    done
done
if [ -e orphans ] ; then
    echo There are $(ls orphans | wc -l) orphaned samples and/or \
        jobs in the orphans directory
    echo a '"rm -Rf orphans/*/ ; rm -R orphans"' will delete them
fi
