@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\snake_eyesss\workflow
call activate snake_eyesss
snakemake --config user=jcsimon --profile profiles/config_Mac -s build_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\imports\Jacob\20241101_3
snakemake --config user=jcsimon --profile profiles/config_Mac -s build_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\imports\Jacob\20241101_4
snakemake --config user=jcsimon --profile profiles/config_Mac -s build_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\imports\Jacob\20241025_1
snakemake --config user=jcsimon --profile profiles/config_Mac -s build_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\imports\Jacob\20241025_2
snakemake --config user=jcsimon --profile profiles/config_Mac -s build_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\imports\Jacob\20241025_3
PAUSE