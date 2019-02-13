#lang racket/base
#| UNDER CONSTRUCTION. The intermediate N space is unneeded. U space is now W space, C space is now Z space, and
   the mapping is performed in one step, w = J(z) + K(z).

   The functionality of this file will mostly reside in transform.rkt. This file probably will keep the record structure
   stuff, and that will be its purpose, factoring out aerodynamics.

Describes an airfoil as a vector of complex numbers in the user's coordinate system (U space),
   a 'normalized' system that ensures the chord is parallel to the x axis, and centered on the origin (N space),
   and a 'cylinder' space, containing the cylinder to which the user's airfoil will be mapped (C space). |#
(provide(all-defined-out))
(require racket/math
         racket/vector
         "circles.rkt"
         "joukowsky.rkt"
         "powell.rkt")
(module+ test
  (require plot 
           "../airfoils.rkt"
           "test-private.rkt"))

;; Our goal, our output, is the geometry struct
(struct geometry
 ; User input and airfoil characteristics in U space
 (U ; vector of complex defining airfoil in user's coordinates, CCW from trailing edge
  nU ; number of points in U (and in vectors to be defined elsewhere in structure)
  iLE ; index of leading edge
  chord-line ; complex # from LE to TE
  c/4 ; complex, quarter chord point
  chord ; real, magnitude of chord line
  (user-chord #:mutable) ; By default, projection of chord line onto x axis, aka (real-part chord-line)
  alpha-0L ; zero lift angle of attack
 
  ; Normalized coordinates and parameters to transform to them from U space
  N
  translate ; translation from U space to N space
  rotate ; rotation from user x axis to normal x axis
  
  ; Circularized coordinates & transform parameters
  C ; vector of complex on circle centered about c0
  c0; center of circle onto which airfoil is mapped.;
  r ; radius of circle onto which airfoil is mapped.
  a ; real, approximately 1/4 of chord
  dT/dc ; derivative of transform, which will be used to calculate velocity
  
  ; Geometric quantities that are calculated, if needed, and cached
  [S #:mutable])) ; distance CCW from trailing edge - don't use mutablility, it's an optimization
                  ; to defer calculation of S, and avoid multiple calculations

;; Our input is the user's coordinates
(define(make-geometry U)
 (let*[(nU(vector-length U))
       (iLE(findLE U))
       (chord-line(-(vector-ref U 0)
                    (vector-ref U iLE)))
       (c/4 (- (vector-ref U 0)
               (* 3/4 chord-line)))
       (chord(magnitude chord-line))
       (user-chord(real-part chord-line))
       (translate(-(vector-ref U 0)(/ chord-line 2)))
       (rotate(-(angle chord-line)))
       (N(for/vector[(u U)]
           (u->n rotate translate u)))
       (a(/ chord 4.0))
       (c0(powell(λ(c)(airfoil-fit-error N a c iLE nU))
                 -0.02+0.02i 0.01 0.001))
       (r(radius a c0))
       (C(make-vector nU))
       (K(make-vector nU))
       (dT/dc(make-vector nU))]
;   ; Adjust a so chord of fitted, user airfoils match
;   (set! a(* a (/ chord(magnitude(-(J(vector-ref C 0)
;                                   (vector-ref N iLE))))))
;   ; Refine c0 using new value of a.
;   (set! c0(powell(λ(c)(airfoil-fit-error N a c iLE nU))
;                  c0 0.005 0.0005))
   ; Fill in C & K
   (vector-set! C 0 a)
   (vector-set! K 0 0)
   (let loop[(i 1)]
     (if(< i iLE)
        (begin
          (vector-set! C i (c-closest-to-n-upper c0 r a (vector-ref N i)))
          (loop(add1 i)))
        (let loop[(i i)]
          (cond[(< i nU)
                (vector-set! C i(c-closest-to-n-lower c0 r a (vector-ref N i)))
                (vector-set! K i(-(vector-ref N i)
                                  (J a(vector-ref C i))))
                (loop(add1 i))]))))
   ; Fill in dT/dc
   (vector-set! dT/dc 0 0.)
   (let loop[(i 1)]
     (if(< i(sub1 nU))
        (begin
          (vector-set! dT/dc i(+(dJ/dc a(vector-ref C i))
                                (dk/dc (vector-ref K (add1 i))
                                       (vector-ref K i)
                                       (vector-ref K(sub1 i))
                                       (vector-ref C(add1 i))
                                       (vector-ref C i)
                                       (vector-ref C(sub1 i)))))
          (loop(add1 i)))
        (begin
          (vector-set! dT/dc i 0.)
          (geometry U nU iLE chord-line c/4 chord user-chord (+ rotate(angle(- a c0))) N translate rotate C c0 r a dT/dc #f))))))


;; Cast a user point in normal coordinates
(define(u->n rotate translate u)
  (*(- u translate) ; Translate first
    (exp(* 0+i rotate)))); positive rotation is CCW
(module+ test
  ; Imagine a ray from 0+3i to 4. If u->n works correctly, transformed ray should go from -2.5 to +2.5.
  (let[(r(-(angle 4-3i)))
       (t 2+3/2i)]
    (check-complex-=(u->n r t 4)2.5 epsilon)))



#| We need to find a point on a circle in C space that is closest to a point n.
     n   -> R(n) -> |c|=a -> |c-c0|=r
   R(n)  -> |c|=a via flat plate transormation
  |c|=a  -> |c-c0|=r via projecting ray from (origin, c0) through a circle onto r circle
                             Instinct says c0 ^^^^^^^^^^ good place for experimentation later |#
(define (c-closest-to-n-upper c0 r a n)
  (let*[(theta(acos(/(real-part n)2.0 a)))
        (ca(* a(exp(* 0+i theta))))]
    (+ c0(* r(exp(* 0+i(angle(- ca c0))))))))
(define (c-closest-to-n-lower c0 r a n)
  (let*[(theta(- (* 2 pi)(acos(/(real-part n)2.0 a)))); only difference with prior function
        (ca(* a(exp(* 0+i theta))))]
    (+ c0(* r(exp(* 0+i(angle(- ca c0))))))))
(module+ test
  ; TE maps to TE
  (check-complex-= 0.25(c-closest-to-n-upper -0.01 0.26 0.25 0.5)epsilon)
  (check-complex-= 0.25(c-closest-to-n-lower -0.01 0.26 0.25 0.5)epsilon)
  ; LE maps to LE
  (check-complex-= -0.27(c-closest-to-n-upper -0.01 0.26 0.25 -0.5)epsilon)
  (check-complex-= -0.27(c-closest-to-n-lower -0.01 0.26 0.25 -0.5)epsilon)
  ; Midpoint maps to near +/- ai
  (check-complex-= 0+0.25i(c-closest-to-n-upper -0.01 0.26 0.25 0)0.01)
  (check-complex-= 0-0.25i(c-closest-to-n-lower -0.01 0.26 0.25 0)0.01))
           
;; We need a measure of how well our Joukowski airfoil fits the user's airfoil
(define(airfoil-fit-error N a c0 iLE nU)
  (let[(r(radius a c0))]
    (let loop[(err 0.0)
              (i 1)]; when things are working right, error at i == 0, nU-1 is zero
      (if (< i iLE)
          (let[(n(vector-ref N i))]
            (loop(+ err(sqr(magnitude(-(vector-ref N i)
                                       (J a(c-closest-to-n-upper c0 r a n))))))
                 (add1 i)))
          (let loop[(err err)
                    (i i)]
            (if (< i nU)
                (let[(n(vector-ref N i))]
                  (loop(+ err(sqr(magnitude(-(vector-ref N i)
                                             (J a(c-closest-to-n-lower c0 r a n))))))
                       (add1 i)))
                err))))))         
(module+ test
  ; Our test airfoil is Joukowski airfoil centered on c0=-.01. Our fitted J-foil is centered on c0=0. 
  (let*[(n 21)
        (N1 (joukowsky 0.25 -0.01 n))    ; Airfoil being "fitted" in normal coordinates
        (N2 (joukowsky 0.25  -0.0001 n))]; Nearly a flat plate - true flat plate caused problems at LE
    ; Our procedure should be able to find better points than merely assuming that theta1 == theta2
    (check-true(> (for/fold
                      [(err 0.0)]
                      [(n1 N1)
                       (n2 N2)]
                    (+ err (sqr(magnitude(- n1 n2)))))
                  (airfoil-fit-error N1 0.25 -0.0001(findLE N1)n))))
    ; Feel good test - does error surface look right?
    (display(plot3d (surface3d (λ(x y)
                                 (airfoil-fit-error (vector-map (λ(n)(- n 1/2))naca-4412); normalize the naca
                                                    0.25 (make-rectangular x y)(findLE naca-4412)(vector-length naca-4412)))) 
                    #:x-min -0.04
                    #:x-max  0.00
                    #:y-min  0.00
                    #:y-max  0.04))
    (printf "~nYou should see a reasonable 3D plot of the airfoil fit function.~n"))

;; The leading edge is estimated by finding the point in U with the minimum real value
(define(findLE U)
  (let loop[(i 1)
            (last-x(real-part(vector-ref U 0)))]
    (let[(this-x(real-part(vector-ref U i)))]
               (if(> this-x last-x)
                  (sub1 i) ; Last point had minimum x value, and thus is on leading edge
                  (loop(add1 i)this-x)))))

;; Estimate derivative of complex function @cm (the middle c)
(define (dk/dc kl km kr cl cm cr)
  (let[(dkl(- kl km))
       (dkr(- kr km))
       (dcl(- cl cm))
       (dcr(- cr cm))]
    (/(-(* dkl(sqr dcr))(* dkr(sqr dcl)))
      (-(* dcl(sqr dcr))(* dcr(sqr dcl))))))
(module+ test
  (let*[(c1(exp 0+3.0i)); A series of points on a circle, similar to what we'll need in practice
        (c2(exp 0+3.1i))
        (c3(exp 0+3.2i))
        (f1(λ(z)(+ 1. (* 2. z)))); Linear functions should be easy
        (f2(λ(z)(-(* 2.(sqr z))(* 3. z)4.))); So should quadratics
        (f3(λ(z)(+(* 3.14 z z z)(* 2.71 z z)(* 1.41 z)0+i))); Cubics?
        (f4(λ(z)(/ 1.0 z z))); As in first term of Laurent series approach
        (f5 sin)]; And some weird shit
    (check-complex-= 2.
                     (dk/dc(f1 c1)(f1 c2)(f1 c3)c1 c2 c3)epsilon)
    (check-complex-=(-(* 4. c2)3.)
                    (dk/dc(f2 c1)(f2 c2)(f2 c3)c1 c2 c3)epsilon)
    (check-complex-=(+(* 3. 3.14 c2 c2)(* 2. 2.71 c2)1.41)
                    (dk/dc(f3 c1)(f3 c2)(f3 c3)c1 c2 c3)0.04); Good fit ~ 0.4% error
    (check-complex-=(/ -2. c2 c2 c2)
                    (dk/dc(f4 c1)(f4 c2)(f4 c3)c1 c2 c3)0.06); Fair fit ~ 3% error
    (check-complex-=(cos c2)
                    (dk/dc(f5 c1)(f5 c2)(f5 c3)c1 c2 c3)0.001))); Good fit ~ 0.2% error

;; Calculate distance from TE, if required
(define (get-S geo) ; overrides default, implements caching behaviour
  (let[(S (geometry-S geo))]
    (if S S ; we've asked for S before
        (let[(N (geometry-N geo))]
          (define old-n (vector-ref N 0))
          (define s 0.)
          (set-geometry-S! geo 
                           (for/vector [(n N)]
                             (set! s (+ s (magnitude (- n old-n))))
                             (set! old-n n)
                             s))
          (geometry-S geo)))))
 
(module+ test
  (let[(my-geo (make-geometry #(4 0+3i 4)))]
    (check-true (not (geometry-S my-geo)))
    (check-true (eqv? (get-S my-geo)
                      #(0 5 10))) ; exercises first time code
    (check-true (eqv? (get-S my-geo)
                      #(0 5 10))))) ; exercises cache code

;; Unit test code and exports
(module+ test
  #| Since this whole unit is based on the Joukowski airfoil, it ought to be able to analyze one.
     Our estimate of a will introduce error, and our estimate of dT/dc must be able to
     compensate for it |#
  (let*[(a 0.25)
        (e -0.01+0.01i)
        (n 21)
        (r(radius a e))
        (test-foil(joukowsky a e n))
        (foil-analysis(make-geometry test-foil))
        (C(JAGC a e n))]
    (printf"Should see two or three errors -- LE accuracy needs work, pretty good elsewhere.~n")
    (for[(dt/dc(geometry-dT/dc foil-analysis))
         (c C)]
      (check-complex-= dt/dc (dJ/dc a c) 0.08))))

(module+ main
  ; Running this module from DrRacket will run all tests.
  (require (submod ".." test)))

