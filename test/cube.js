"use strict"

var tape = require("tape")
var uniq = require("uniq")
var ch = require("../ich")

function createHypercube(dimension, numSubdivs) {
  var initSamples = new Array(numSubdivs+2)
  initSamples[0] = 0.0
  initSamples[numSubdivs+1] = 1.0
  for(var i=1; i<=numSubdivs; ++i) {
    var theta = i / (numSubdivs+2)
    initSamples[i] = theta
  }
  var points = initSamples.map(function(x) {
    return [x]
  })
  for(var i=1; i<dimension; ++i) {
    var npoints = []
    for(var j=0; j<points.length; ++j) {
      for(var k=0; k<initSamples.length; ++k) {
        var p = points[j].slice()
        p.push(initSamples[k])
        npoints.push(p)
      }
    }
    points = npoints
  }
  return npoints
}

function validateCube(t, dimension) {
  var verts = createHypercube(dimension, 3)

  //Shuffle vertices

  var x = new Array(dimension)
  for(var i=0; i<dimension; ++i) {
    x[i] = 0
  }
  verts.unshift(x)
  for(var i=0; i<dimension; ++i) {
    var y = x.slice()
    y[i] = 1
    verts.unshift(y)
  }

  var hull = ch(verts, true)

  var hullverts = []
  for(var i=0; i<hull.length; ++i) {
    t.equals(hull[i].length, dimension, "check facet: " + hull[i].join())
    hullverts.push.apply(hullverts, hull[i])
  }
  uniq(hullverts)


  for(var i=0; i<hullverts.length; ++i) {
    var p = verts[hullverts[i]]
    var onboundary = false
    for(var j=0; j<dimension; ++j) {
      if(p[j] === 0 || p[j] === 1) {
        onboundary = true
      }
    }
    t.ok(onboundary, "point on boundary: " + p.join())
  }
}

tape("cube-2d", function(t) {
  for(var i=0; i<10; ++i) {
    validateCube(t, 2)
  }
  t.end()
})

tape("cube-3d", function(t) {
  for(var i=0; i<10; ++i) {
    validateCube(t, 3)
  }
  t.end()
})

tape("cube-4d", function(t) {
  validateCube(t, 4)
  t.end()
})