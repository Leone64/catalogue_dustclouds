# A small catalogue of dust clouds within 1.25 kpc of the Sun
---

A repository to store data and scripts relevant to my [Bachelor Thesis](https://leone64.github.io/bachthesis/).

---

### FuncDef

A module which defines the commonly used functions in Mass calculation and plotting.

---

### ErrorCalc

A script which calculates the masses and the uncertainties of all the clouds defined in ```cloud-data.csv```.

Note: this current version was improved by [Simon Waldmann](https://github.com/waldini1) to enable it to run on 6 cores simultaneously. My PC only has 6 cores, and crashes on start of the script.

---

### CSV-Files

```cloud-data.csv``` stores the parameters of the Clouds in the Catalogue, as presented in Tab. 3.1 of the [Thesis](https://leone64.github.io/bachthesis/Bachelor_Thesis.pdf).

```results.csv``` is the output of the ```ErrorCalc.py``` script.

---

### Plotting_Script

An example plotting script. Cloud 18 of the catalogue is plotted in 4 different ways, similar to Figs. 3.1 and 3.3 of the Thesis, see the interactive version [here](https://leone64.github.io/bachthesis/dropdown/).

The output of the plotting script is stored in ```savedplots/```.