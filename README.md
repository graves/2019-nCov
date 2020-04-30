# 2019-nCov
SIR model implemented in Smalltalk for modelling COVID-19 transmission.

Check out [my blog](https://graves.github.io/covid-19-SIR-model.html) for the how and why this was built.

```Smalltalk
| model view |
model := SIR initializeWithDurationOfRecovery: 7.5 populationSize: 8399000 r0: 2.6.
model := model solveODEWithNinfected: 1 nRecovered: 0 nDays: 180.
view := model chartFromDate: 'March 1, 2020'.
model exportViewAsPNG.
model exportViewAsHTML.
```

![Example PNG](https://i.imgur.com/6V1lpix.png)
