# Group 23 MATLAB code
`Nov 25, 2026 · @AV`
### How to run
Two scripts do everything; run them from the inner `Group-23-main\Group-23-main` folder.
- `Simple`: analyses the Canvas example data, then runs the theoretical fuel analysis and the load study; all plots land in one tabbed window (several minutes, 72 load points).
- `RunExperimentAnalysis`: analyses every run in Data/ExperimentRuns.csv, calibrates the model and predicts injection timings; works without Simple.
- Results go to new folders each run: `FuelAnalysisResults/, LoadStudyResults/, ExperimentResults/<date_time>/`.

### Measured-data analysis (Simple.m)

Turns a measured pressure trace into work, efficiency and combustion timing.
- **`Load`**: pressure file (crank angle, pressure in bar), 3600 points per cycle (0.2°), 69 cycles.
- **`processMeasuredPressure`**: pegs each cycle to 1.00 bar at intake BDC (−180°), averages the cycles, smooths with a 21-point cubic Savitzky–Golay filter.
- **`CylinderVolume`**: slider-crank volume V(θ); swept volume 0.722 L, compression ratio 21.5.
- **`Work and IMEP`**: area of the p–V loop per cycle (trapz), then power and efficiency.
- **`Heat release`**: aROHR per cycle with dynamic γ from calcGamma (NASA cp/cv, gas changes from air to exhaust), averaged; CA10/50/90 from the normalised cumulative heat between −20° and 80°.
- **`addIdealDieselCycle`**: ideal Diesel cycle drawn on the pV plots as a reference (γ = 1.35, cut-off ratio 2).
- **`KPIs`**: use placeholder flows (fuel 0.3 g/s); replace with lab values.
### Theoretical fuel model
Simulates the cylinder from equations, with no measured pressure, to compare fuels and timings.
- **`faFuels`**: 6 fuels: Diesel (= B7), FT diesel, OME3, HVO, B7-OME3 12.5 % and 25 % (blends built from B7 and OME3 with mixing rules).
- **`faSimulate`**: closed cycle from −180° to 180°; energy equation with NASA gas properties; heat released along a prescribed burn curve after an ignition delay of 0.8 ms × 52.2/CN. No heat loss or gas exchange, so efficiencies (about 62–68 %) are idealised.
- **`faRunAnalysis`**: 6 fuels × SOI 14, 16, 18, 21° bTDC at equal fuel energy (1032 J/cycle) → KPI table and plot tabs.
- **`faRunLoadStudy`**: solves the fuel mass (fzero) that gives 1.5, 3.0 and 4.5 bar gross IMEP → 72 operating points.
### Experimental results analysis
Processes any number of lab runs, each with its own fuel, load, SOI and flows.
- `Data/ExperimentRuns.csv`: one row per run: data file, fuel, engine, load, SOI, intake pressure and temperature, fuel and air flow, CO2/CO/HC/NOx, brake power. Empty cells give NaN; Include = 0 skips a run.
- `faAnalyzeExperiment`: same steps as Simple.m, but with that run's SOI, intake conditions and fuel. Outputs net and gross IMEP, cycle-to-cycle variation (COV), peak pressure, start of combustion, ignition delay, CA10/50/90, efficiency, ISFC, CO2 per kWh (carbon balance), λ and specific NOx/CO/HC; brake values when brake power is entered.
- `faPlotExperiment`: tabs for pressure, pV, heat release, γ and temperature, KPIs against SOI, and the KPI table.
### Model validation and timing prediction
Tunes the theoretical model on one measured timing, then predicts the others.
- `faValidateModel`: groups runs by fuel, engine and load; calibrates faSimulate on the run nearest 18° bTDC using its ignition delay, burn shape (CA10/50/90), heat-release fraction and compression pressure. faSimulate itself is unchanged.
- `Prediction`: 14–21° bTDC at constant gross IMEP (fuel mass re-solved per SOI, as in the lab); the best timing is the highest efficiency within optional pressure limits.
- `Validation`: other measured timings of the same group are compared with the predictions (CA50, peak pressure, fuel mass, efficiency).
- `faPlotValidation`: calibration plot, prediction plots and tables.
### Limitations and tips
- Only the example data is analysed so far; its fuel flow (0.3 g/s), intake temperature (293 K) and pegging pressure (1 bar) are placeholders.
- The "bsfc/bsCO2/bsNOx" in Simple.m use p–V work, so they are indicated values, not brake values.
- Theoretical efficiencies are idealised: compare trends, not absolute values. HVO and FT diesel are nearly identical in the model (same formula).
- Before adding lab files, check the format matches (crank angle in column 1, pressure in bar in column 2, cycles starting at −360°); otherwise change the settings at the top of RunExperimentAnalysis.m.
- Every run creates a new results folder; delete old ones now and then.
