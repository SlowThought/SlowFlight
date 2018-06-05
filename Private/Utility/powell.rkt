#lang racket/base
(require math/number-theory
         racket/math)
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

;; Find minimum of real function of vector of reals of arbitrary dimension.
;; Initially, search is along normal unit vectors. Search results are used to
;; modify the search vectors to approximate a gradient descent.
(define (powell f x0 h epsilon max)
  ; f(x) is function to be minimized - x is a vector
  ; x0 is starting position, a vector
  ; h is initial step size along search vector
  ; epsilon is step size limit, if we move less than epsilon in a given direction, we stop
  ; max is the max number of iterations of the algorithm
  (define n (vector-length x0))
  (define n+1 (add1 n))
  (define maxi (* max n))
  ; We use xs and s as ring buffers
  (define (index-s i)(with-modulus n (mod i)))
  (define (index-x i)(with-modulus n+1 (mod i)))
  ; Build initial set of search vectors
  (define s (build-vector n (λ (i)
                              (define v (make-vector n 0.))
                              (vector-set! v i 1.)
                              v)))
  ; History of search
  (define xs (make-vector n+1 x0))
  (define fs (make-vector n+1 (f x0)))
  
  ; Search along initial directions
  (vector-set! fs 0 (f x0))
  (for [(i (in-range n))]
    (define-values (x* f*)(explore-along-s f (vector-ref xs i)(vector-ref s i) h epsilon))
    (vector-set! xs (add1 i) x*)
    (vector-set! fs (add1 i) f*))
  ; Adjust directions based on previous search results (excluding last one, as that direction
  ; is fruitless, for the moment).
  (for [(i (in-range n maxi))]
    (vector-set! s (index-s i)(vector-map - (vector-ref xs (index-x (sub1 i)))
                                            (vector-ref xs (index-x (- i n)))))
    (define-values (x* f*)(explore-along-s f (vector-ref xs (index-x i))(vector-ref s (index-s i)) h epsilon))
    (vector-set! xs (index-x (add1 i)) x*)
    (vector-set! fs (index-x (add1 i)) f*))
  (values (vector-ref xs (index-x maxi))
          (vector-ref fs (index-x maxi))))
;; powell helper functions
(require racket/vector)
(define (explore-along-s f xm s h epsilon)
  (let*[(xl (vector-map (λ (x s)(- x (* s h))) xm s))
        (xr (vector-map (λ (x s)(+ x (* s h))) xm s))
        (fl (f xl))
        (fm (f xm))
        (fr (f xr))
        (ds (find-minimum fl fm fr
                        (- h)0.  h))
        (x* (vector-map (λ (x s) (+ x (* s ds))) xm s))
        (f* (f x*))]
    (values x* f*)))
              
(module+ test
  ; Our test function is a hopefully well behaved quartic, with square  terms that will skew things enough to keep
  ; the search interesting
  (define (f xs)
    (let[(v (vector-ref xs 0))
         (w (vector-ref xs 1))
         (x (vector-ref xs 2))
         (y (vector-ref xs 3))
         (z (vector-ref xs 4))]
      (+ (* 3 (- v 1)(- v 1)(- v 1)(- v 1))
         (* 2 (- w 2)(- w 2)(- w 2)(- w 2))
         (* 4 (- x 3)(- x 3)(- x 3)(- x 3))
         (* 7 (- y 4)(- y 4)(- y 4)(- y 4))
         (* 5 (- z 6)(- z 6)(- z 6)(- z 6))
         (* v x)
         (* w z))))
  ; Check values derived from testing involving lots of printf statements which have since been edited out. Consider
  ; adding debug argument to optionally enable noisy output
  (define-values (x* f*)(powell f #(0 0 0 0 0) 0.1 0.01 2))
  (check-equal? x* #(0.7186233703322971 1.975163313924521 3.384259279958702 4.775251092524473 5.056951459999327))
  (check-= f* 19.00948681490754 epsilon))
(module+ main
  (require(submod".."test)))




  
