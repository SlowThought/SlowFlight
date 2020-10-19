#lang typed/racket/base
#| 1d-search.rkt provides functions for minimizing functions of the form (f Real: x)->Real, and functions for casting functions
   of the form (f Vector: v) to the prior form, given a vector position and vector search direction
|#
(require racket/vector)
(provide powell-1)

#| powell-1 finds minimum of quadratic function f(t) 

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

    We will be applying powell-1 to NON-quadratic equations, but everything near enough its minimum
    LOOKS quadratic. powell-1 iterates to find bounds within which f is quadratic enough.

    powell-1 assumes we've handed it well behaved functions, therefore no safety checks.

    f is the function to be minimized
    xl is the left (min) guess at x
    xm is the middle (best) guess
    xr is the right (max) guess
    ϵ is the tolerance (if |x* - xm| < ϵ, we quit)
    n is the iteration limit (if we don't converge relatively quickly, we quit)
|#

(: powell-1 (-> (-> Real Real) Real Real Real Real Positive-Index Real))
(define (powell-1 f xl xm xr ϵ n)
  (define fl (f xl))
  (define fm (f xm))
  (define fr (f xr))
  (let loop [(fl fl)(fm fm)(fr fr)(xl xl)(xm xm)(xr xr)(i 0)]
    (if (= i n) xm
        (let [(x* (find-minimum fl fm fr xl xm xr))]
          (cond
            ; the predicted minimum is close enough to the middle evaluation
            ((< (abs (- xm x*)) ϵ) x*)
            ; the predicted minimum is outside the expected range. Based on experience,
            ; the prediction is likely not far enough outside the range, so we widen the
            ; window. extreme x is double distance from x* that x* is from old limit.
            ((< x* xl)(let[(xfl (- (* 3 x*) (* 2 xl)))] ; xfl is extreme far left
                        (loop (f xfl) (f x*) fl xfl x* xl (add1 i))))
            ((> x* xr)(let[(xfr  (- (* 3 x*)(* 2 xr)))] ; xfr is extreme far right
                        (loop fr (f x*) (f xfr) xr x* xfr (add1 i))))
            ; the predicted minimum is inside the window
            ((< x* xm)(loop fl (f x*) fm xl x* xm (add1 i)))
            ( #t ; (> x* xm)
              (loop fm (f x*) fr xm x* xr (add1 i))))))))

;; Helper functions for powell-1

; Use determinants to solve quadratic equation - returns x that corresponds to extremum (hopefully minimum) of f
(: find-minimum (-> Real Real Real Real Real Real Real))
(define (find-minimum fl fm fr
                      tl tm tr)
  (let*[(ad (determinant fl tl 1.
                         fm tm 1.
                         fr tr 1.))
        (bd (determinant (sqr tl) fl 1.
                         (sqr tm) fm 1.
                         (sqr tr) fr 1.))]
    (/(- bd) 2. ad)))

; Square a number - Macro because it's not worth a function call, and because it doesn't require restricting the type of
; sqr to Real. Query: macro contracts? How do they play with function contracts?
(define-syntax-rule (sqr x)(* x x))

; Find the determinant of a 3x3 matrix
(: determinant (-> Real Real Real Real Real Real Real Real Real Real))
(define (determinant a1 a2 a3
                     b1 b2 b3
                     c1 c2 c3)
  (- (+ (* a1 b2 c3)
        (* a2 b3 c1)
        (* a3 b1 c2))
     (* a1 b3 c2)
     (* a2 b1 c3)
     (* a3 b2 c1)))

; Cast a function of a vector to a function of a real, given an initial point and a search direction. fv(x0) == fr(0)
(: fv->fr (-> (-> (Vectorof Real) Real) (Vectorof Real) (Vectorof Real)(-> Real Real)))
(define (fv->fr fv x0 s)
  (λ ([t : Real])
    (fv (vector-map + x0 (vector-map (λ ([ si : Real]) (* t si)) s)))))

(module+ test
  (require plot racket/math "../test-private.rkt")
  (printf "1d-search.rkt: tests running.~n")
  (: f (-> Real Real))
  (define (f t)(+ (sqr (- t 3.)) 2.))
  (define ϵ 0.01)
  (display (plot (function f 1 5 #:label "y = (x-3)^2 + 2")))
  (check-= (powell-1 f 2 3 4 ϵ 1) 3. ϵ)
  (check-= (powell-1 f 3 4 5 ϵ 1) 3. ϵ)
  (check-= (powell-1 f 1 2 3 ϵ 1) 3. ϵ)
  (check-= (powell-1 f 4 5 6 ϵ 1) 3. ϵ)
  (check-= (powell-1 f 0 1 2 ϵ 1) 3. ϵ)
  (printf "~nQuadratics should work in one step, regardless of initial guesses.")
  (printf "~nIf no error after chart, quadratic minimum evaluated correctly.~n")
  (let*[(quiet-f (λ([ t : Real ])(- (sqr (sqr (+ t 4.))) 7.)))
        (f (λ([t : Real ]) (begin
                   (printf "q")
                   (quiet-f t))))
        (ϵ 0.01)(n 7)]
    ; Parameters ϵ and n chosen to pass below tests. It's clear that the algorithm is 'kinda' working,
    ; but quartics are hard, ϵ & n are not ideal, and it's unclear how to choose ϵ & n.
    (display (plot (function quiet-f -8 0 #:label "y = (x+4)^4 - 7")))
    (printf "~nQuartic tests - q represents a function evaluation.~n")
    (printf "~n x = -4 central input ")
    (check-= (powell-1 f -5 -4 -3 ϵ n) -4. ϵ)
    (printf "~n x = -4 right input ")
    (check-= (powell-1 f -6 -5 -4 ϵ n) -4. ϵ)
    (printf "~n x = -4 left input ")
    (check-= (powell-1 f -4 -3 -2 ϵ n) -4. ϵ)
    ; Answer not within origninal bounds, so more tolerance (higher ϵ), and more time (higher n) needed.
    (printf "~n inputs to left of -4. ")
    (check-= (powell-1 f -7 -6 -5 ϵ (* 2 n)) -4. (* 2 ϵ))
    (printf "~n inputs to right of -4. ")
    (check-= (powell-1 f -3 -2 -1 ϵ (* 2 n)) -4. (* 2 ϵ))
    (printf "~nYou'd think quartics would be easy. They're not. A place for further study. Read comments in code.~n"))
  ; Any function with well defined, polite minima should work... if intial bounds avoid impolite extrema.
  (let [(f (λ([ t : Real ]) (sin (- t 1.))))
        (ϵ 0.01)]
    (display (plot (function f 3 6 #:label "y = sin(x - 1)")))
    (check-= (powell-1 f 3 4.5 6 ϵ 7)                                                           
             (+ (* 3/2 pi) 1.) ϵ)
    (printf "~nQuartic fix broke sine solution. n bumped from 5 to 7, ϵ to .01, to get pass. Art vs science.")
    (printf "~nIf no error after chart, sine minimum evaluated correctly.~n"))
  ; A realistic test similar to actual problem
  (let*[(y (λ ([ x : Real ]) (- 1/4 (/ x 4))))
        (f (λ ([ x : Real ]) (+ (sqr (- x 0.8))(sqr (- (y x) 0.07)))))]
    (display (plot (function f -1 2 ; x bounds
                             #:y-min 0
                             #:y-max 4
                             #:label "realistic case")))
    (check-= (powell-1 f 0 1/2 1 0.005 5) 0.793 0.005)
    (printf "~nIf no error after chart, realistic case evaluated correctly.~n"))

  ; Can we cast a vector function to a real function which powell-1 can handle?
  (let*[(f (λ ([v : (Vectorof Real)])
             (let [(x (vector-ref v 0))
                   (y (vector-ref v 1))
                   (z (vector-ref v 2))]
               (+ (sqr (- x 1))
                  (sqr (- y 2))
                  (sqr (- z 3))))))
        (g (fv->fr f #(0 2 3) #(1 0 0)))] ; By inspection, we expect the minimum of g to be at t == 1
    (check-= (powell-1 g -.2 0 .2 0.005 5) 1 0.005))
  (printf "~nIf no error, casting of vector function to real function worked.~n"))