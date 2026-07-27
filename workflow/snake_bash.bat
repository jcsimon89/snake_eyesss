@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\snake_eyesss\workflow
call activate snake_eyesss

timeout /t 1 /nobreak

snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS261\fly_009 --unlock
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS261\fly_009




PAUSE

