@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\snake_brainsss\workflow
call activate snake_brainsss
snakemake --config user=jcsimon --profile profiles/config_Mac -s build_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\brainsss_imports\20240531
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\brainsss\JS68_x_JS117\fly_005
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\brainsss\JS68_x_JS117\fly_006
PAUSE