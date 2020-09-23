# Regression


[![Build Status](https://travis-ci.org/shaiq681/Regression.svg?branch=master)](https://travis-ci.org/shaiq681/Regression)

## Introduction

An R Package to implement ls() style regression.

## Installation

```r
devtools::install_github("shaiq681/Regression")
```

## Usage

- For creating a new Regression Object

```r
linreg_obj <- linreg$new(formula= Sepal.Length~Petal.Width, data = iris, use.qr = F)
```



