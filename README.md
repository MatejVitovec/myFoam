# myFoam

## Sestavení

Nejdřív načtěte prostředí kompatibilní instalace OpenFOAM (cestu upravte podle své instalace), potom spusťte skript z kořene projektu:

```bash
source /cesta/k/OpenFOAM/etc/bashrc
./Allwmake
```

Skript postupně sestaví knihovny v `src/`, testovací a pomocné programy a aktuální řešitele v `applications/solvers/`. Při chybě sestavení skončí hned a vypíše neúspěšný příkaz. Archivované řešitele v `OLDSOLVERS/` záměrně nesestavuje.
