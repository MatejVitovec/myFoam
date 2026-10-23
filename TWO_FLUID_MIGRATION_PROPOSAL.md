# Navrh migrace two-fluid solveru

Tento navrh je pouze ke kontrole. Zatim nejsou zmeneny zadne zdrojove soubory solveru.

## Cile migrace

Referencni implementaci predstavuje [rk3TwoFluidTurbulentFoam](applications/solvers/multiphase/rk3TwoFluidTurbulentFoam/rk3TwoFluidTurbulentFoam.C): pouziva `CompressibleTwoFluidFoam::compressibleTwoFluid`, jeho konzervativni stav a fluxy, novy drag model a `CompressibleTwoFluidFoam::fvc::div`.

Ostatni solvery stale pouzivaji puvodni `TwoFluidFoam::twoFluid`:

| Solver | Zvlastni chovani, ktere musi zustat zachovano |
| --- | --- |
| `dbFoam` | Explicitni RK3 + drag |
| `rk3PrimitiveTwoFluidFoam` | Aktualizace primitivnich promennych pomoci `calculateDeltaPrimitive.H` |
| `lusgsTwoFluidFoam` | LUSGS sweep a konservativni rezidua |
| `lusgsTwoFluidTurbulentFoam` | LUSGS + turbulence + virtual-mass cleny |
| `lusgsTwoFluidBoundedFoam` | LUSGS sweep s vlastni logikou pro ohraniceni alpha |
| `lusgsDbFoam` | LUSGS + kondenzacni zdrojove cleny |
| `wetSteamDbFoam` | RK3 + kondenzace/nukleace/rust kapek |
| `lusgsTwoFluidWetSteamFoam` | LUSGS + kondenzacni zdroje a latentni teplo |

## Spolecne zmeny

Ve zdrojovem `.C` a `createFields.H` nahradit stare typy a hlavicky API novymi podle referencniho solveru. Typicky jde o:

```diff
-#include "twoFluid.H"
-#include "twoFluidConvectiveFlux.H"
-#include "twoFluidFvc.H"
+#include "compressibleTwoFluid.H"
+#include "compressibleTwoFluidConvectiveFlux.H"
+#include "compressibleTwoFluidFvc.H"

-TwoFluidFoam::twoFluid fluid(mesh);
-twoFluidConservative& conservative = fluid.conservative();
+CompressibleTwoFluidFoam::compressibleTwoFluid fluid(mesh);
+CompressibleTwoFluidFoam::compressibleTwoFluidConservative& conservative = fluid.conservative();
```

Dale je potreba pouzit novy typ `compressibleTwoFluidConvectiveFlux`, prevadet `TwoFluidFoam::fvc::div` na `CompressibleTwoFluidFoam::fvc::div`, a zachovat poradi aktualizaci stavu jako v referenci (`correctBoundaryCondition`, `correctThermo`, `correctInterfacialPressure`, `correctConservative`). Konkretni vyrazy fluxu a zdroju se musi porovnat po jednotlivych solverech; nejde bezpecne provest globalni textovou nahradou.

V `Make/options` nahradit include cestu `src/twoFluid/lnInclude` cestou `src/phaseSystemModels/compressibleTwoFluid/lnInclude` a knihovnu `-lmyTwoFluid` knihovnou `-lmyCompressibleTwoFluid`. Turbulentni varianty navic potrebuji include cestu a knihovnu turbulence modelu podle Make souboru referencniho solveru. Ostatni existujici zavislosti (napr. `myCondensationModels`, `myBoundFoam`, turbulence a thermophysical knihovny) ponechat, pokud je solver pouziva.

## Zmeny specificke pro skupiny

1. `dbFoam`: nejblizsi prvni kandidat pro migraci zakladniho explicitniho konzervativniho reseni. Zachovat jeho RK3 koeficienty, drag a ukladane veliciny; porovnat zejmena tlakove fluxy fazi s referencnim solverem.
2. `rk3PrimitiveTwoFluidFoam`: rozhodnout, zda zachovat primitivni aktualizaci a pouze napojit nove `fluid`/flux API, nebo prevest solver na konzervativni aktualizaci. Druha varianta meni numerickou formulaci, proto ji nepovazuji za automatickou soucast migrace.
3. `lusgsTwoFluidFoam`, `lusgsTwoFluidTurbulentFoam`, `lusgsTwoFluidBoundedFoam`: upravit rezidua a sweep helpery spolu s hlavnim `.C`; zachovat LUSGS matici, relaxaci, virtual mass a omezeni alpha. Zmeny se tykaji take `forwardSweep*.H`, `backwardSweep*.H` a pripadnych dalsich helperu.
4. `lusgsDbFoam`, `wetSteamDbFoam`, `lusgsTwoFluidWetSteamFoam`: ponechat kondenzacni/nukleacni zdroje mimo jadro fluidu a zachovat jejich bilance hmoty, hybnosti a energie. Po zmene fluxu znovu zkontrolovat znamenka zdroju a vazbu prumeru kapek na drag model.

## Doporuceny postup po schvaleni

1. Migrovat `dbFoam` a zbuildit ho jako prvni referencni port bez zmeny schematu.
2. Samostatne migrovat `rk3PrimitiveTwoFluidFoam` po potvrzeni, ze ma zustat primitive solver.
3. Migrovat kazdy LUSGS solver vcetne jeho sweep helperu a otestovat sestaveni.
4. Nakonec migrovat tri wet-steam/condensation varianty a overit zachovani zdrojovych clenu na existujicich tutorialech.

Pro kazdou skupinu porovnat nejmene startovni sestaveni, beh kratkeho tutorialu a bilance hmoty/energie. V teto fazi jde o navrh rozsahu, nikoli o potvrzeni numericke ekvivalence.