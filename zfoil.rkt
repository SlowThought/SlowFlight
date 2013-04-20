#lang racket/base
;; Provide Racket subset most useful to engineers
(require racket/math
         racket/vector)
(provide(all-from-out racket/base
                      racket/math
                      racket/vector))

;; Provide contracts for zfoil functions to be exported
(require"airfoils.rkt"
        "private/geometry.rkt"
        "private/inviscid.rkt"
        "private/plots.rkt"
        racket/contract)
(provide[contract-out
         ; From airfoils.rkt
         (read-dat-file(-> string? (vectorof complex?)))
         (naca-4412 (vectorof complex?))
         (naca-66_2-415 (vectorof complex?))
         ; From geometry.rkt
         (geometry-alpha-0L(-> geometry? real?))
         (geometry-chord(-> geometry? real?))
         (geometry-U(-> geometry? (vectorof complex?)))
         (make-geometry(-> vector? ; of complex
                           geometry?))
         ; From inviscid.rkt
         (V(-> real? geometry? vector?)) ; of complex
         (MagV(-> vector? geometry? vector?))
         (CL(->* (real? geometry?)
                 ([vectorof real?])
                 real?))  
         (Cp(-> real? real?))]
         ; From plots.rkt (I have to figure out what 'renderer-tree' is in contract-speak)
         plot-foil
         plot-vs-x)

       
         