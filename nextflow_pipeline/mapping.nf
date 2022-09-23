#!/usr/bin/env nextflow

Channel.from( ['A', 1, 2], ['B', 3, 4] )
    .map{ it -> (it[0], it[1]+it[2]) }
    .view()

Channel.from( ['A', 1, 2], ['B', 3, 4] )
    .map{ it -> [ it[1], it[2] ] }
    .view()

