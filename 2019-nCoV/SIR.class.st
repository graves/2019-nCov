Class {
	#name : #SIR,
	#superclass : #Object,
	#instVars : [
		'dt',
		'durationOfRecovery',
		'gamma',
		'populationSize',
		'r0',
		'beta',
		'odeValues',
		'view'
	],
	#category : #'2019-nCoV'
}

{ #category : #initialization }
SIR class >> initializeWithDurationOfRecovery: nDays  populationSize: nPersons r0: aFloat [
	"comment stating purpose of class-side message"
	"scope: class-variables  &  class-instance-variables"	
			
	| sir |
	
	sir := self basicNew.
	sir dt: 1.0.
	sir durationOfRecovery: nDays.
	sir populationSize: nPersons.
	sir r0: aFloat.
	
	^ sir
]

{ #category : #accessing }
SIR >> beta [
	beta ifNil: [
			self beta: self r0 * self gamma / self populationSize ].
	^ beta
]

{ #category : #accessing }
SIR >> beta: aFloat [
	beta := aFloat
]

{ #category : #'as yet unclassified' }
SIR >> chartFromDate: aStartDateString [
	"comment stating purpose of instance-side message"

	"scope: class-variables  &  instance-variables"

	| startDate b ds1 ds2 ds3 |
	startDate := aStartDateString asDate.
	b := RTGrapher new.
	ds1 := RTData new.
	ds1 label: 'Susceptible'.
	ds1 noDot.
	ds1 points: self susceptible.
	ds1 connectColor: Color blue.
	ds1 y: [ :v | v y ].
	ds1 x: [ :v | v x ].
	ds1 interaction toggleDataset.
	ds1 interaction popup
		text: [ :v | 
			v key y asInteger asString , ' susceptible on '
				, (startDate + v key x day) asDate asString ].
	b add: ds1.
	ds2 := RTData new.
	ds2 label: 'Infected'.
	ds2 noDot.
	ds2 points: self infected.
	ds2 connectColor: Color red.
	ds2 y: [ :v | v y ].
	ds2 x: [ :v | v x ].
	ds2 interaction toggleDataset.
	ds2 interaction popup
		text: [ :v | 
			v key y asInteger asString , ' infected on '
				, (startDate + v key x day) asDate asString ].
	b add: ds2.
	ds3 := RTData new.
	ds3 label: 'Recovered'.
	ds3 noDot.
	ds3 points: self recovered.
	ds3 connectColor: Color green.
	ds3 y: [ :v | v y ].
	ds3 x: [ :v | v x ].
	ds3 interaction toggleDataset.
	ds3 interaction popup
		text: [ :v | 
			v key y asInteger asString , ' recovered on  '
				, (startDate + v key x day) asDate asString ].
	b add: ds3.
	b addDecorator: (RTCursorFollower new color: Color gray).
	b axisX
		title: '';
		labelRotation: -30;
		labelConversion: [ :v | (startDate + v day) asDate ].
	b axisY
		title: 'Population';
		noDecimal.
	b legend right.
	view := b build view.
	^ view
]

{ #category : #accessing }
SIR >> dt [
	^ dt
]

{ #category : #accessing }
SIR >> dt: aFloat [
	dt := aFloat
]

{ #category : #accessing }
SIR >> durationOfRecovery [
	^ durationOfRecovery
]

{ #category : #accessing }
SIR >> durationOfRecovery: nDays [
	durationOfRecovery := nDays
]

{ #category : #'as yet unclassified' }
SIR >> exportViewAsHTML [
	"comment stating purpose of instance-side message"
	"scope: class-variables  &  instance-variables"	
			
	RTHTML5Exporter new export: view
]

{ #category : #'as yet unclassified' }
SIR >> exportViewAsPNG [
	"comment stating purpose of instance-side message"
	"scope: class-variables  &  instance-variables"	
			
	view canvas buildMorph extent: 1000@500; exportAsPNG
]

{ #category : #accessing }
SIR >> gamma [
	gamma ifNil: [ 
		self gamma: 1 / durationOfRecovery ].
	^ gamma
]

{ #category : #accessing }
SIR >> gamma: aFloat [
	gamma := aFloat
]

{ #category : #'as yet unclassified' }
SIR >> infected [
	"comment stating purpose of instance-side message"
	"scope: class-variables  &  instance-variables"	
			
	^ self odeValues collectWithIndex: [ :each :idx | Point x: idx y: (each at: 2) ]

]

{ #category : #initialization }
SIR >> initialize [
	"comment stating purpose of instance-side message"
	"scope: class-variables  &  instance-variables"	
			
	dt := 1.0.
]

{ #category : #accessing }
SIR >> odeValues [
	^ odeValues
]

{ #category : #accessing }
SIR >> odeValues: anArray [
	odeValues := anArray
]

{ #category : #accessing }
SIR >> populationSize [
	^ populationSize
]

{ #category : #accessing }
SIR >> populationSize: nPersons [
	populationSize := nPersons
]

{ #category : #accessing }
SIR >> r0 [
	^ r0
]

{ #category : #accessing }
SIR >> r0: aFloat [
	r0 := aFloat
]

{ #category : #'as yet unclassified' }
SIR >> recovered [
	"comment stating purpose of instance-side message"
	"scope: class-variables  &  instance-variables"	
			
	^ self odeValues collectWithIndex: [ :each :idx | Point x: idx y: (each at: 3) ]

]

{ #category : #accessing }
SIR >> solveODEWithNinfected: nInfected nRecovered: nRecovered nDays: nDays [

	| days system stepper solver state |
	
	days := nDays asFloat.
	
	self gamma.
	self beta.
	
	system := PMExplicitSystem block: [ :x :t| |c|
		c := Array new: 3.
		c at: 1 put: (self beta negated) * (x at: 1) * (x at: 2).
		c at: 2 put: (self beta * (x at: 1) * (x at: 2)) - (self gamma * (x at: 2)).
		c at: 3 put: self gamma * (x at: 2).
		c
		].
	
	stepper := PMRungeKuttaStepper onSystem: system.
	solver := (PMExplicitSolver new) stepper: stepper; system: system; dt: dt.
	state := { populationSize . nInfected . nRecovered }.
	self odeValues: ((0.0 to: days by: dt) collect: [ :t| state := stepper doStep: state
																			time: t stepSize: dt ]).
	^ self
]

{ #category : #'as yet unclassified' }
SIR >> susceptible [
	"comment stating purpose of instance-side message"
	"scope: class-variables  &  instance-variables"	
			
	^ self odeValues collectWithIndex: [ :each :idx | Point x: idx y: (each at: 1) ]

]

{ #category : #accessing }
SIR >> view [
	^ view
]

{ #category : #accessing }
SIR >> view: anObject [
	view := anObject
]
