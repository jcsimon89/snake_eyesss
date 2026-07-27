@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\snake_eyesss\workflow
call activate snake_eyesss

timeout /t 1 /nobreak

snakemake --config user=jcsimon --profile profiles/config_Mac_buildfly -s build_fly.smk --directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\imports\Jacob\20260716


PAUSE