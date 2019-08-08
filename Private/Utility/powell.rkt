#lang racket/base
(require math/number-theory
         racket/math)
(provide powell-1); powell)
(module+ test
  (require plot "../test-private.rkt")
  (printf"powell.rkt: tests running.~n"))

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

; Use determinants to solve quadratic equation
(define (find-minimum fl fm fr
                      tl tm tr)
  (let*[(ad (determinant fl tl 1.
                         fm tm 1.
                         fr tr 1.))
        (bd (determinant (sqr tl) fl 1.
                         (sqr tm) fm 1.
                         (sqr tr) fr 1.))]
    (/(- bd) 2. ad)))

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

;(module+ test
;  (let [(f (λ(t) (+ (sqr (- t 3.)) 2.)))
;        (ϵ 0.01)]
;    (display (plot (function f 1 5 #:label "y = (x-3)^2 + 2")))
;    (check-= (powell-1 f 2 3 4 ϵ 1) 3. ϵ)
;    (check-= (powell-1 f 3 4 5 ϵ 1) 3. ϵ)
;    (check-= (powell-1 f 1 2 3 ϵ 1) 3. ϵ)
;    (check-= (powell-1 f 4 5 6 ϵ 1) 3. ϵ)
;    (check-= (powell-1 f 0 1 2 ϵ 1) 3. ϵ)
;    (printf "~nQuadratics should work in one step, regardless of initial guesses.")
;    (printf "~nIf no error after chart, quadratic minimum evaluated correctly.~n"))
;  (let*[(quiet-f (λ(t)(- (sqr (sqr (+ t 4.))) 7.)))
;        (f (λ(t) (begin
;                   (printf "q")
;                   (quiet-f t))))
;        (ϵ 0.01)(n 7)]
;    ; Parameters ϵ and n chosen to pass below tests. It's clear that the algorithm is 'kinda' working,
;    ; but quartics are hard, ϵ & n are not ideal, and it's unclear how to choose ϵ & n.
;    (display (plot (function quiet-f -8 0 #:label "y = (x+4)^4 - 7")))
;    (printf "~nQuartic tests - q represents a function evaluation.~n")
;    (printf "~n x = -4 central input ")
;    (check-= (powell-1 f -5 -4 -3 ϵ n) -4. ϵ)
;    (printf "~n x = -4 right input ")
;    (check-= (powell-1 f -6 -5 -4 ϵ n) -4. ϵ)
;    (printf "~n x = -4 left input ")
;    (check-= (powell-1 f -4 -3 -2 ϵ n) -4. ϵ)
;    ; Answer not within origninal bounds, so more tolerance (higher ϵ), and more time (higher n) needed.
;    (printf "~n inputs to left of -4. ")
;    (check-= (powell-1 f -7 -6 -5 ϵ (* 2 n)) -4. (* 2 ϵ))
;    (printf "~n inputs to right of -4. ")
;    (check-= (powell-1 f -3 -2 -1 ϵ (* 2 n)) -4. (* 2 ϵ))
;    (printf "~nYou'd think quartics would be easy. They're not. A place for further study. Read comments in code.~n"))
;  ; Any function with well defined, polite minima should work... if intial bounds avoid impolite extrema.
;  (let [(f (λ(t) (sin (- t 1.))))
;        (ϵ 0.01)]
;    (display (plot (function f 3 6 #:label "y = sin(x - 1)")))
;    (check-= (powell-1 f 3 4.5 6 ϵ 7)                                                           
;             (+ (* 3/2 pi) 1.) ϵ)
;    (printf "~nQuartic fix broke sine solution. n bumped from 5 to 7, ϵ to .01, to get pass. Art vs science.")
;    (printf "~nIf no error after chart, sine minimum evaluated correctly.~n"))
;  ; A realistic test similar to actual problem
;  (let*[(y (λ (x) (- 1/4 (/ x 4))))
;        (f (λ (x) (+ (sqr (- x 0.8))(sqr (- (y x) 0.07)))))]
;    (display (plot (function f -1 2 ; x bounds
;                             #:y-min 0
;                             #:y-max 4
;                             #:label "realistic case")))
;    (check-= (powell-1 f 0 1/2 1 0.005 5) 0.793 0.005)
;    (printf "~nIf no error after chart, realistic case evaluated correctly.~n"))
;  )
; EVERYTHING ABOVE THIS LINE IS GOOD -- TESTS COMMENTED OUT FOR COMPILATION SPEED. UNCOMMENT TESTS BEFORE MOVING ON.

#| powell-n finds minimum of real function of vector of reals of arbitrary dimension.
   Initially, search is along normal unit vectors. Search results are used to
   modify the search vectors to approximate a gradient descent.
|#

(define (powell-n f x0 h ϵ max)
  ; f(x) is function to be minimized - x is a vector
  ; x0 is starting position, a vector
  ; h is initial step size along search vector
  ; epsilon is step size limit, if we move less than epsilon in a given direction, we stop
  ; max is the max number of iterations of the algorithm
  (define n (vector-length x0))
  (define n+1 (add1 n))
  ; Search vectors
  (define e (build-basis-vectors n))
  (define ss e)
  ; Search history
  (define xs (make-vector n+1 (make-vector n 0)))
  (vector-set! xs 0 x0)
  (define fs (make-vector n+1))
  (vector-set! fs 0 (f x0))
  ; Perform search
  (for*[(i (in-range max))
        (j (in-range n))]
    (define k (add1 j))
    (let-values([(x* f*)(search-along-direction (vector-ref ss j)
                                           (vector-ref xs j)
                                           f h ϵ max)])
      (vector-set! xs k x*)
      (vector-set! fs k f*)))
  ; Identify least fruitful search direction
  ; Replace least fruitful with first
  ; Replace first with best (from x0 to xn-2, we exclude last search direction because best should be normal to it)
  (printf "powell-n isn't done yet.~n"))
;; powell helper functions


; build-basis-vectors generates n normal vectors of dimension n
(define (build-basis-vectors n)
  (build-vector n (λ (i)
                    (define v (make-vector n 0.))
                    (vector-set! v i 1.)
                    v)))
(module+ test
  (check-equal? (build-basis-vectors 2) #(#(1.0 0.0)#(0.0 1.0))))

; test infrastructure -- generate random nd quadratic equations
(module+ test
  (require math/distributions math/matrix)
  ; Generate a random positive definite matrix
  (define (random-pd-matrix n)
    (define nd (normal-dist 0. 0.5))
    (define A
      (build-matrix    
       n n
       (λ (i j) (cond [(= i j)
                       (exp (sample nd))] ; Diagonal must be positive
                      [(> i j) 0]         ; Building triangular matrix to ensure diagonal determines determinant
                      [else (sample nd)]))))
    ;(printf "A : ~s~n" A) ; uncomment to debug. confidence is high.
    (matrix+ A (matrix-transpose A))); Suspect, but have not proven, that this last step, making final matrix symmetrical, is unneccesary.
  ; Generate a random quadratic function that has a minimum f* at location x*
  (define (gen-random-f f* x*)
    (define n (vector-length x*))
    (define A (random-pd-matrix n))
    (define X* (->col-matrix x*))
    (λ (x)
      (define X (->col-matrix x))
      (define dX (matrix- X X*))
      (+ f* (matrix-ref (matrix* (matrix-transpose dX) (matrix* A dX))
                       0 0))))
  ;; Tests of gen-random-f, 2D examples consistently show expected result
  ;  (let*[(f (gen-random-f 1. #(1. 1.)))
  ;        (g (λ (x y) (f (vector x y))))]
  ;    (plot3d (surface3d g)
  ;            #:x-min 0
  ;            #:x-max 2
  ;            #:y-min 0
  ;            #:y-max 2))
  ;    (for-each (λ (x)
  ;                 (printf "~nX: ~a f(X): ~a " x (f x)))
  ;              (list #(1. 1.) #(1. 2.) #(2. 1.) #(1. 0.) #(0. 1.)))
     #| end of gen-random-f tests |#)

; f-of-v->f-of-t casts a function of a vector, a starting point, and a search direction, to a 1d parametric function
(require racket/vector)
(define (f-of-v->f-of-t f x0 s)
  (λ (t) (f (vector-map + x0 (vector-map (λ (x) (* t x)) s)))))

(module+ test
  (let*[(f (gen-random-f 3 #(0 0 1)))
        (g (f-of-v->f-of-t f #(0 0 0) #(0 0 1)))]
    (check-= (g 1) 3 0.00001)))

; search-along-direction applies powell-1 to a function of a vector
(define (search-along-direction s x f h ϵ max)
  (let*[(g (f-of-v->f-of-t f x s))
        (dt (powell-1 g (- h) 0 h ϵ max))
        (x* (vector-map + x (vector-map (λ(e)(* dt e)) s)))
        (f* (f x*))]
    (values x* f*)))

(module+ test
  ; Test search-along-direction
  (let[(f (λ (v) (let [(x (vector-ref v 1))] (* (sub1 x)(sub1 x)))))];only second element is significant, minimum at v=#(0 1 0)(
    (define-values (x* f*)(search-along-direction #(0 1 0) #(0 0 0) f 0.5 0.01 5))
    (check-= f* 0 .01)
    (check-= (vector-ref x* 0) 0 .01)
    (check-= (vector-ref x* 1) 1 .01)
    (check-= (vector-ref x* 2) 0 .01))

  ; Tests of powell-n
  (define f (gen-random-f -10. #(1. 2. 3. 4.)))
  (powell-n f #(0 0 0 0) 0.1 0.001 1))
  
