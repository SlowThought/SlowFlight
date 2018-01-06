#lang racket/base
(require racket/math)
(provide find-minimum); powell)
(module+ test
  (require plot "../test-private.rkt")
  (printf"powell.rkt: tests running.~n"))

#| Find minimum of quadratic function f(t) 

   f(t) = at^2 + bt +c

   Minimum (ok, extreme, may be min or max, but we'll be building well behaved functions) when
   derivative is zero.

   f' = 2at + b = 0. tmin = -b/2a.

   We are solving the linear system

     | tl^2 tl 1||a| |fl|
     | tm^2 tm 1||b}=|fm|
     | tr^2 tr 1||c| |fr|

    We don't need c for our minimum calculation. We use method of determinants to find a and b.
    We don't actcually calculate a or b, as both involve dividing by the same determinant d, and
    so the ds cancel in the final calculation.
  |#
(define (find-minimum fl fm fr
                      tl tm tr)
  (let*[(ad (determinant fl tl 1.
                         fm tm 1.
                         fr tr 1.))
        (bd (determinant (sqr tl) fl 1.
                         (sqr tm) fm 1.
                         (sqr tr) fr 1.))]
    (/(- bd) 2. ad)))
;; Helper functions for find-minimum
; Square a number
(define (sqr x)(* x x))
; Find the determinant of a 3x3 matrix
(define (determinant a1 a2 a3
                     b1 b2 b3
                     c1 c2 c3)
  (- (+ (* a1 b2 c3)
        (* a2 b3 c1)
        (* a3 b1 c2))
     (* a1 b3 c2)
     (* a2 b1 c3)
     (* a3 b2 c1)))
(module+ test
  (let [(f (λ(t) (+ (sqr (- t 3.)) 2.)))]
    (display (plot (function f 1 5 #:label "y = (x-3)^2 + 2")))
    (check-= (find-minimum (f 2.)(f 3.)(f 4.)
                              2.    3.    4.)
             3. epsilon)
    (printf "~nIf no error after chart, quadratic minimum evaluated correctly.~n"))
  (let [(f (λ(t) (- (sqr (sqr (+ t 4.))) 7.)))]
    (display (plot (function f -8 0 #:label "y = (x+4)^4 - 7")))
    (check-= (find-minimum (f -5.)(f -3.)(f -1.)
                              -5.    -3.    -1.)
             -4. epsilon)
    (printf "~nIf no error after chart, quartic minimum evaluated correctly.~n"))
  ; Any function with well defined, polite minima should work... iteration should let us get close enough.
  (let [(f (λ(t) (sin (- t 1.))))]
    (display (plot (function f 3 6 #:label "y = sin(x - 1)")))
    (check-= (let loop [(i 5)(h 1.)(tm 5.)]
               (if (= i 0) tm
                   (let [(tl (- tm h))
                         (tr (+ tm h))]
                    (loop (sub1 i)(/ h 2.)(find-minimum (f tl)(f tm)(f tr)
                                                           tl    tm    tr)))))                                                           
             (+ (* 3/2 pi) 1.) 0.00001)
    (printf "~nIf no error after chart, sine minimum evaluated correctly.~n")))

;;; Find minimum of real function of a complex variable
;(define (powell f z0 h epsilon)
;  ; First we move in h direction
;  (let[(z1(+ z0(* h(find-minimum(f(- z0 h))
;                                (f z0)
;                                (f(+ z0 h))))))]
;    ; Then we move normal to it
;    (let*[(i(* 0+i h))
;          (z2(+ z1(* i(find-minimum(f(- z1 i))
;                                   (f z1)
;                                   (f(+ z1 i))))))
;          ; Then we find a new step that sort of is near the steepest gradient
;          (new-h(- z2 z0))]
;      (if (< (magnitude new-h) epsilon); we're done
;          z2; and we return the point corresponding to the minimum
;          (powell f z2 new-h epsilon))))); else we search some more in a hopefully better direction
;(module+ test
;  (check-= (magnitude(powell (λ(z)(test-function(real-part z)(imag-part z)))-10.0+100.0i 1.0 epsilon))
;           0.0 epsilon))
(module+ main
  (require(submod".."test)))




  
