@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\snake_eyesss\workflow
call activate snake_eyesss

snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS259\fly_001 --unlock
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS259\fly_001
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS259\fly_002 --unlock
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS259\fly_002
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS259\fly_003 --unlock
snakemake --config user=jcsimon --profile profiles/config_Mac -s preprocess_fly_dev.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS140_x_JS259\fly_003

PAUSE
