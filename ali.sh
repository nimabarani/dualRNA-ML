

for log_path in /home/nima/projects/def-lpenacas/nima/newDual/scratch/raw_reads/*/*/*/star/star_Log.final.out
do
    dirname="$(dirname "$log_path")" # get the parent directory of the file path
    srid="$(basename "$(dirname "$dirname")")"
    type="$(basename "$(dirname "$(dirname "$dirname")")")"
    experiment_name="$(basename "$(dirname "$(dirname "$(dirname "$dirname")")")")"

    cp $log_path /home/nima/projects/def-lpenacas/nima/newDual/scratch/logs/$experiment_name/$type/$srid.out
done