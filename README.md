# Dockerized version of PC+ RISM model

## Setting up

```
git clone git@github.com:MTS-Strathclyde/pcplus-docker.git
cd pcplus-docker
docker build . -t rism-worker
docker run -p 3002:3002 rism-worker
```

Check that it works (in another terminal):
```
> curl localhost:3002/health
OK
```


## Running calculations

Compute properties of ethane in water:

```
curl --request POST \
  --url http://localhost:3002/calculate \
  --header 'content-type: application/json' \
  --data '{
	"SMILES": "CC",
	"solvent": "water",
	"closure": "pse3",
	"tolerance": 1.0e-4
}'
```

expected result:
```
{
  "converged": true,
  "results": {
    "PC_exchem": -0.492866397439224,
    "PC_plus_exchem": 1.2247003246756893,
    "PMV": 86.98072471046082,
    "closure_exchem": 11.94775020401471,
    "solute_solvent_potential_energy": -4.581540443508729
  }
}
```

All energy results are in `kcal/mol` units, partial molar volume is expresse in cubic Angstroms.

More details: https://github.com/MTS-Strathclyde/PC_plus


## API

URL: `/calculate`

METHOD: `POST`

PAYLOAD:
```
{
	"SMILES": <Valid molecule in SMILES format>
	"solvent": <Any option out of the list bellow>,
	"closure": <kh|pse2|pse3|hnc>,
	"tolerance": <Positive float, values ~1.0e-5 recommended>
}
```

SUPPORTED SOLVENTS:
```
['acetonitrile', 'benzene', 'bromobenzene',
'carbtet', 'chloroform', 'CS2', 'water',
'cyclohexane', 'decane', 'dichloroethane',
'diethylether', 'dmso', 'ethylacetate', 'heptane',
'isooctane', '0.1M NaCl solution', '0.2M NaCl
solution', '0.3M NaCl solution', '0.4M NaCl
solution', '0.5M NaCl solution', 'octanol',
'olive_oil', 'toluene', 'xylene']
```




