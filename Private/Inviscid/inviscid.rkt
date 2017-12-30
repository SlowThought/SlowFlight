#lang racket/base
;;; Calculate circulation about a cylinder, resulting velocity about an airfoil, Cl, Cm
(provide(all-defined-out))
(require racket/math
         "geometry.rkt"
         "quad.rkt"
         "zvector.rkt")

(module+ test
  (require  plot
            racket/vector
            "circles.rkt"
            "test-private.rkt"
            "../airfoils.rkt")
  (printf"inviscid.rkt running tests.~n")
  ; generate a Joukowski airfoil, for which we have an analytic solution to compare against
  (define j-pts 25)
  (define j-thetas (JAGC 1. .01 j-pts)) ; Joukowski Airfoil Generating Circle
  (define j-foil (joukowsky 1. .01 j-pts))
  (define j-geo (make-geometry j-foil))
  (define (j-V theta alpha)
    #| from http://slow-flight.blogspot.com/2013/04/how-complex-numbers-make-things-simple.html
       theta corresponds to a point on the JAGC |#
    (let*[(a 1.)
          (e -.01)
          (r (radius a e))
          (w (-(* r(exp(* 0+i theta)))e))
          (dPhi/dw (-(exp(* 0-1i alpha))
                     (/(*(sqr r)(exp(* 0+1i alpha)))
                       (sqr(+ w e)))
                     (/(* 0-2.02i(sin alpha))
                       (+ w e))))]
      ; V = dPhi/dw / df/dw
      (/ dPhi/dw (- 1 (/(sqr w)))))))
          

;; Given angle of attack, geometry, return velocity distribution as a vector
;; We build a complex potential function superimposing free stream, doublet & vortex
(define(V alpha geo)
  (let*[; Extract geometry
        (alpha-C(- alpha (geometry-rotate geo))); alpha in "cylinder space"
        (r(geometry-r geo)); this turns out to be cooefficient for doublet term
        (a(geometry-a geo))
        (c0(geometry-c0 geo))
        (last(sub1(geometry-nU geo)))
        ; Free stream
        (v-inf(exp(* 0+i alpha-C)))
        ; Doublet
        (k1(* -1 r r v-inf))
        ; Calculate circulation term
        (k2(-(/(* r r v-inf)(- a c0))
             (*(conjugate v-inf)(- a c0))))
        ; Build potential function
        (dPhi/dc(λ(c)(+(conjugate v-inf)(/ k1(sqr(- c c0)))(/ k2(- c c0)))))]
    (let[(ret(for/vector
                 [(c(geometry-C geo))
                  (dt/dc(geometry-dT/dc geo))]
               ; Use chain rule to find velocity in user space
               (conjugate(/(dPhi/dc c)dt/dc))))]
      ; Trailing edge is almost certainly NAN - we force it to be a stagnation point
      (vector-set! ret 0 0.+0.i)
      (vector-set! ret last 0.+0.i)
      ret)))
      
; Translate velocity vectors to velocity magnitudes, including sign
(define(MagV vs geo)
  (define zs(geometry-U geo))
  (define n(vector-length zs))
  (define ms(make-vector n))
  (for[(i(in-range 1(sub1 n)))]
    (let[(ds(-(vector-ref zs (add1 i))
              (vector-ref zs (sub1 i))))]
      (vector-set! ms i(/(zdot ds(vector-ref vs i))
                         (magnitude ds)))))
  ms)

(module+ test
  (let*[(iLE (geometry-iLE j-geo))
        (Vs (V (geometry-alpha-0L j-geo) j-geo))
        (ms (MagV Vs j-geo))
        (velocity-tolerance .15)
        (magnitude-tolerance .1)]
    ; easy case, 0 angle of attack
    (for[(i(in-range 1 (sub1 j-pts)))
         (v Vs)
         (m ms)]
      (cond
        [(< i iLE)
         (check-complex-= (vector-ref Vs i) 1. velocity-tolerance)
         (check-= (vector-ref ms i) -1. magnitude-tolerance)]
        ; we expect LE to be screwed up on a flat plate
        [(> i iLE)
         (check-complex-= (vector-ref Vs i) 1. velocity-tolerance)
         (check-= (vector-ref ms i) 1. magnitude-tolerance)]))
    ; now we step it up.
    (define alpha 0.1)
    (set! Vs (V alpha j-geo))
    (set! ms (MagV Vs j-geo))
    (for[(i(in-range 0 j-pts))
         (x j-foil)
         (v Vs)
         (m ms)]
      (cond[(or(= i 0)
               (= i iLE)
               (= i (sub1 j-pts)))
            ; We know these are exceptional points, and do not expect good answers
            (void)]
           [(< i iLE)
            (define theta (acos(/(real-part x)2.)))
            (define true-v (/(+(sin(- theta 0.1))(sin 0.1))(sin theta)))
            (check-complex-= v true-v velocity-tolerance)
            (check-= (- m) true-v magnitude-tolerance)]
           [(> i iLE)
            (define theta (-(* 2. pi)(acos(/(real-part x)2.))))
            (define true-v (/(+(sin(- theta 0.1))(sin 0.1))(sin theta))) (check-complex-= v true-v velocity-tolerance)
            (check-= (abs m) true-v magnitude-tolerance)]))))   
    
; Coefficient of Pressure, for incompressible flow, V_inf == 1
(define(Cp v)
  (- 1 (sqr v)))

;; Find lift coefficient
(define (CL alpha geo 
            [vs (MagV (V alpha geo) geo)])           
  (define zs(geometry-U geo))
  (define n (sub1(vector-length zs)))
  (let loop[(i 0)(CX 0.)(CY 0.)]
    (if(< i n)
       (let*[(z1(vector-ref zs i))
             (x1(real-part z1))
             (y1(imag-part z1))
             (z2(vector-ref zs(add1 i)))
             (x2(real-part z2))
             (y2(imag-part z2))
             (v1(vector-ref vs i))
             (v2(vector-ref vs(add1 i)))]
         (loop(add1 i)
              (+ CX (integrate v1 v2 x1 x2 2))
              (+ CY (integrate v1 v2 y1 y2 2))))
       (-(/ (+ (* CX (cos alpha))
               (* CY (sin alpha)))
            (geometry-user-chord geo))))))

(module+ test
  (define lift-tolerance 0.1)
  ; We first test the integration routine by providing a bogus velocity distribution
  ; THIS CHECK ASSUMED VALUE OF J-POINTS. MUST WRITE MORE GENERIC VERSION
  ;  (check-= (CL 0. j-geo #(-1. -1.1 -1.2 -1.3 -1.4 -1.5 -1.6 -1.7 -1.8 -1.9 -1.991 0. 1.
  ;                              1.  1.   1.   1.   1.   1.   1.   1.   1.   1.)) 2/3 lift-tolerance)
  ; We test a real airfoil w/ free stream || zero lift line
  (let-values[((zs comment)(read-dat-file"../Dat/r1046.dat"))]
    (let*[(geo(make-geometry zs))
          (a0L(geometry-alpha-0L geo))
          (vs (MagV(V a0L geo)geo))]
      (display (plot (points (vector-map 
                              (λ(z v)(vector (real-part z)(-(Cp v))))
                              zs vs))
                     #:title (format "Alpha = ~a" a0L)
                     #:y-label "-Cp"))
      (check-= (CL a0L geo vs) 0. lift-tolerance)
      ; Now we check at positive lift
      (define alpha1 (+ a0L 0.1))
      (set! vs (MagV (V alpha1 geo) geo))
      (display (plot (points (vector-map 
                              (λ(z v)(vector (real-part z)(-(Cp v))))
                              zs vs))
               #:title (format "Alpha = ~a" alpha1)
               #:y-label "-Cp"))
      (check-= (CL alpha1 geo vs)(* 0.2 pi) lift-tolerance))))

;; Find moment coefficient
; z0 is the point about which we're calculating the moment.
; z1 and z2 define the interval over which we're integrating.
; v1, v2 are velocities at z1, z2
(define (dM z0 z1 z2 v1 v2)
  (let*[(delta-z(- z2 z1)); surface coordinate axis
        (|s|(magnitude delta-z))
        (n^ (/ delta-z 0-i |s|)); unit vector normal to surface
        ; assume v is of form v=a+bs between z1, z2, s goes from 0 to |s|
        (a v1)
        (b (/(- v2 v1)|s|))
        ; z is of same form on line segment between z1, z2, z=c+ds
        (c z1)
        (d (/ delta-z |s|))
        ; dM = (z-z0) x n^ p ds
        ; we assume rho=2, p0=0, so Bernoulli pressure relation becomes p=-v^2,
        ; a function of s, as is z
        ; dM = -[(c + ds - z0) x n^](a + bs)^2 ds
        ;    = -[(c - z0) x n^ + ds x n^](a^2 + 2abs + b^2s^2)ds
        (e (zcross (- c z0) n^))
        (f (zcross d n^))
        (g (sqr a))
        (h (* 2. a b))
        (i (sqr b))]
    ; dM = -(e + fs)(g + hs + is^2)ds
    ;    = -[eg + (fg + eh)s + (fh + ei)s^2 + fis^3]ds
    ;  M = -[egs + 1/2(fg + eh)s^2 + 1/3(fh + ei)s^3 + 1/4(fi)s^4]
    (+ (* e g |s|)
       (* (/ 2.)(+ (* f g)(* e h))(sqr |s|))
       (* (/ 3.)(+ (* f h)(* e i))(expt |s| 3))
       (* (/ 4.) f i (expt |s| 4)))))
        
(module+ test
  ; With constant velocity over interval, moment about midpoint should be zero
  (check-=  0.(dM .5 0. 1. 1. 1.) epsilon)
  ; Pressure, as defined, should be a suction. Moment from above distribution,
  ; measured about beginning of interval, should be -.5
  (check-= .5(dM 0. 0. 1. 1. 1.) epsilon))


; We sum total moment and non-dimensionalize
(define (CM alpha geo 
            [vs (MagV (V (- alpha (geometry-alpha-0L geo)) geo) geo)])
  (let*[(zs  (geometry-U geo))
        (n   (sub1(vector-length zs)))
        (c/4 (geometry-c/4 geo))]
    (let loop[(i 0)(M 0.)]
      (if(< i n)
         (loop(add1 i)
              (+ M (dM c/4 
                       (vector-ref zs i)
                       (vector-ref zs (add1 i))
                       (vector-ref vs i)
                       (vector-ref vs (add1 i)))))
         (/ M (geometry-user-chord geo))))))

(module+ test
  ; symmetrical airfoil should have zero pitching moment, measured from
  ; quarter chord, regardless of angle of attack.
  (define moment-tolerance 0.025)
  (check-= 0. (CM 0. j-geo) moment-tolerance)
  (check-= 0. (CM .1 j-geo) moment-tolerance))

         

