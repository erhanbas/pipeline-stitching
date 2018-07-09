# pipeline-stitching
computational pipeline that creates yml file that will be consumed by the render

1. Create a folder and clone the repo
```
cd pipe
clone https://github.com/erhanbas/pipeline-stitching.git 
```
2. run main.m in matlab
```
main(brain,[inputfolder,pipelineoutputfolder])
```
## inputs/outputs
    brain : sample name, e.g. '2018-06-14';
    inputfolder [OPTIONAL]: acqusition folder, e.g. sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s',brain);
    pipelineoutputfolder [OPTIONAL]: outputfolder e.g. sprintf('/nrs/mouselight/cluster/sandbox2/%s',brain)

## sample usage
    main('2018-06-14')
    This will create tilebase.cache.yml file in the '/nrs/mouselight/cluster/classifierOutputs/'2018-06-14' folder 
    Then run the render
   
    
## sample render usage
    copy a sample parameter file next to yml file, e.g. 
    copy /nrs/mouselight/cluster/classifierOutputs/2018-07-02/set_parameters_raw.jl /nrs/mouselight/cluster/classifierOutputs/'2018-06-14/parameters_raw.jl
    
    then modify the paths in parameters_raw.jl as follows
    const user="mluser"
    const input="2018-06-14"
    const output="2018-06-14-raw"
    
    then run
    /groups/mousebrainmicro/mousebrainmicro/Software/barycentric4/src/render/src/render parameters_raw.jl

    with default settings, this will create output render in /nrs/mouselight/Users/mluser/2018-06-14
