@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\snake_eyesss\workflow
call activate snake_eyesss
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS015_x_JS251\fly_001
PAUSE