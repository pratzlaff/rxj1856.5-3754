dir=/data/legs/rpete/flight/rxj1856.5-3754

obsids()
{
    grep '^[0-9]' $dir/obsids | cut -f 1 #| tail -1
}

instrument()
{
    local obsid="$1"
    local f=$(ls "$dir"/data/"$obsid"/tg_reprocess/*_evt2.fits 2>/dev/null)

    if [ -z "$f" ]
    then
        grep "$obsid" "$dir"/obsids | perl -anle 'print $F[2]'
        return
    fi

    echo $(detnam "$f")
}

detnam()
{
    local evt2="$1"
    punlearn dmkeypar
    dmkeypar "$evt2" detnam echo+
}
