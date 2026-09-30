# OrdinalCoClust

## To dos

- Coclustering ordinal data function takes in data, number of row and col cluster, number of random inits (20 by default).
- Pour chaque estimation de modele, des initialisations randoms multiples seront implemented en interne de la fonction, et la meilleure solution sera retounree :
	- les probabilites d’appartenances des lignes et colonnes à chaque cluster ligne et cluster colonne,
	- l’affectation des lignes et colonnes aux clusters ligne et colonne,
	- la valeur du critere ICL,
	- les parametres estimés du model.

package R avec la function avec une aide pour votre fonction, incluant un exemple d’utilisation.

Livrable :
	- un package R (en format .tar.gz que je puisse installer sans soucis sur mon Mac),
	- une vignette pour votre package R illustrant l’utilisation de votre package,
	- un oral au cours duquel je vous demanderai de m’expliquer certaines parties de votre code.


EM, CUB

Gibbs, SEM, LBM


|2		|Coclustering (params, z_ij,...)|
|1		|Random init|
|2		|ICL for selecting number of blocks|
|2		|Test|
|4		|Package|
|2		|Example ordinal data|

## Examples

- [Young People Survey](https://www.kaggle.com/cardot/se-young-people-survey/data)

## References

- ordinalClust package.
- CUB package.