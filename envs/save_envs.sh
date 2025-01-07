for env in "CanSig-R" "CanSig-python" "cansig-benchmark" 
do
    mamba env export --from-history -n $env > without_build/$env.yml
    mamba env export -n $env > with_build/$env.yml
done