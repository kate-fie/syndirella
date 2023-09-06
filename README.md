

LIBRARY ENUMERATION BY SYNTHESIS RESTRICTION FOR CAR
=====

> 👩‍🔬 LIBRARY ENUMERATION BY SYNTHESIS RESTRICTION FOR CAR

A pipeline for library enumeration of follow-up compounds from a fragment screen restricted by the set of reactions performed by the XChem 
chemist assisted reaction (CAR) setup and the availability of building blocks.

## Inputs

* Crystallographic data from a fragment screen (i.e. from Fragalysis)
* A set of follow-up compounds

## Outputs

* [x] .csv of analog series from follow-up compounds
* [ ] Fragmenstein placement

# Reactions
(checked means tested)
* [x] Amidation
* [ ] Sp2-Sp2 Suzuki Coupling
* [ ] Buchwald-Hartwig Amination
  * [x] with halides
  * [ ] with OTf
* [ ] Reductive amination
* [x] Formation of urea from two amines
* [x] Sulfonamide Schotten-Baumann with amine (intermolecular)
* [ ] N-nucleophilic aromatic substitution
* [ ] Amide schotten-baumann

# Search
* [x] Manifold superstructure search for building blocks
* [ ] Arthor smallworld search
* [ ] McPairs functional group replacement search
* [ ] Other Manifold searches

```
python main_v2.py -i /Users/kate_fieseler/PycharmProjects/retrievesynthesizable/NCS1/sept4_elab/4Sept_Dani_additional_routes_1_step_edit.csv -r /Users/kate_fieseler/PycharmProjects/retrievesynthesizable/NCS1/sept4_elab/ -u
```