#lang racket/base
;; plots.rkt - Plotting routines for zfoil
(require plot)
(provide (all-defined-out))
(module+ test
  (require racket/math
           "../test-private.rkt");
  (printf"plots.rkt running tests.~n"))

;; plot airfoils. Accepts one or more vectors of complex.
(define (plot-foil . zs)
  (let*[(lol-of-zs (map vector->list zs))
        (xs (tree-map(λ(z)(real-part z)) lol-of-zs))
        (x-min (tree-apply min xs))
        (x-max (tree-apply max xs))
        ; We ASSUME airfoil is thin, and reasonbly centered on x axis
        (y-max (/ (- x-max x-min)2))
        (y-min (- y-max))]
    (plot (cons (axes)
                (map(λ(lozs)(lines(map complex->vector lozs)
                                  #:x-min x-min
                                  #:x-max x-max
                                  #:y-min y-min
                                  #:y-max y-max))
                    lol-of-zs)))))
(module+ test
  (define crude-foil #(1. 0.3+.12i 0+.04i 0.2 1.)) 
  (display (plot-foil crude-foil))
    (printf "~nYou should see a plot of a crude airfoil.~n"))

;; plot coefficients vs x. voz is a vector of complex representing airfoil shape.
;; loloc is one or more lists of real (typically nondimensional cooefficients of some kind).
(define(plot-vs-x voz . loloc)
  (let*[(loz (vector->list voz))
        (lox(map real-part loz))
        (lov(map complex->vector loz))
        (tree-of-v(map(λ(loc)(map(λ(x c)(vector x c))lox loc))
                      loloc))]
    (plot (list (axes)
                (lines lov)
                (map points tree-of-v)))))
(module+ test
  (display(plot-vs-x crude-foil (list 0. 1.5 2. .75 0.)))
  (printf"~nYou should see a simulated plot of Cp vs x.~n"))

;; zfoil often uses arrays of complex, plot expects lists of vectors - if this were C I could 
;; probably speed with a type-cast
(define (complex->vector z)
  (vector (real-part z)(imag-part z)))
;; plot can take multiple sets of data, as a list of lists
; tree-map returns a tree (not really a tree, a list of lists, a lol) 
; of the same shape as lol, with op applied to each value
(define(tree-map op lol)
  (map (λ(l)(map op l))lol))
(module+ test
  (check-equal? (tree-map sqr '((1 -1)
                             (2 -3)))
                '((1 1)
                  (4 9))))
; tree-apply returns a single (probably scalar?) value to an operation applied across a tree
; (not really a tree, a lol).
(define(tree-apply op lol)
  (apply op (map(λ(l)(apply op l))lol)))
(module+ test
  (check-eq?(tree-apply min '((0 1 2 3)
                              (4 5 6 7)))
            0))