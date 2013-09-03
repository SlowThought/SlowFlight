#lang racket/base
#| corner.rkt computes some properties of viscous corner flow, as at the stagnation point of an airfoil.
   This is based on the similarity solution presented by Frank White in "Viscous Fluid Flow", 1st edition.
   If this code works as intended (if I understand Racket), there will be no runtime overhead to modules
   requiring corner.rkt, and I will have successfully used code to document the reasons behind the value of
   a (multiple of a) constant. |#

(provide (all-defined-out))
;(provide (H0 B nu x)      ; Shape factor near stagnation point
;         (delta*0 B nu))  ; Displacement thickness near stagnation point
;         (theta0  B nu x)); Momentum thickness near ...

(require racket/math racket/list)

#| A solution to NS of the form Psi = Bxf(y) is assumed. y and f are non-dimensionalized by the factor (B/nu)^1/2
   Capital letters represent the resulting non-D variables. When all the calculus and algebra are done, we are left
   with the ODE and BCs:

   F''' + FF' + 1 - (F')^2 = 0
   F(0) = F'(0) = 0
   F' -> 1 as Y approaches infinity
   F''(0) = ???

   We shall have to pick a value for F'', and then integrate F''' to find our solution. 
   To increase accuracy, we use derivatives of F''', which is roughly equivalent to 
   what Runga-Kutta does.
   F'''   = (F')^2 - FF'' - 1
   F''''  = 2F'F'' - F'F'' - FF'''
          = F'F'' - FF'''
   F''''' = (F'')^2 - F'F''' - F'F''' -FF''''
          = (F'')^2 - 2F'F''' - FF'''' 
   ...and that's ridiculous enough. |#

(define (F3 F F1 F2)
  (- (sqr F1)(* F F2) 1.))
(define (F4 F F1 F2 F3)
  (- (* F1 F2)(* F F3)))
(define (F5 F F1 F2 F3 F4)
  (- (sqr F2)(* 2. F1 F3)(* F F4)))

(module+ test
  (require "test-private.rkt")
  (printf "corner.rkt running tests...~n")
  (define e 0.001)
  (check-= (F3 0. 0. 1.) -1. e)
  (check-= (F3 999. 1. 0.) 0. e)
  (check-= (F4 0. 0. 1. 0.) 0. e)
  (check-= (F4 999. 1. 0. 0.) 0. e)
  (check-= (F5 0. 0. 1. 0. 0.) 1. e)
  (check-= (F5 999. 1. 0. 0. 0.) 0. e)) 
    
#| A generic integration routine. 
   We have at least 3 terms of Taylor series readily available
   g(x+h) = g(x) + g'(x)h + 1/2 g''(x)h^2 + 1/6 g'''(x)h^3 + ... |#

(define (integrate h gx . derivatives)
  (let loop[(order 1.)
            (divisor 1.)
            (sum gx)
            (derivatives derivatives)]
    (if (null? derivatives)
        sum
        (loop (add1 order)
              (/ divisor (add1 order))
              (+ sum (* (expt h order)(car derivatives)divisor))
              (cdr derivatives)))))

(module+ test
  (check-= (integrate 1. 0. 1.) 1. e)
  (check-= (integrate 1. 0. 1. 1.) 1.5 e)
  (check-= (integrate 1. 0. 1. 1. 1.)(+ 1.5(/ 6.))e))

#| We integrate all the Fs at once, including higher orders as if
   we were using Runga-Kutta |#
(define (integrate-F h F F1 F2) 
  (let*[(F3(F3 F F1 F2))
        (F4(F4 F F1 F2 F3))
        (F5(F5 F F1 F2 F3 F4))]
    (values (integrate h F F1 F2 F3 F4 F5)
            (integrate h F1 F2 F3 F4 F5)
            (integrate h F2 F3 F4 F5))))

#| We need a concept of "big enough" for infinity and "small enough" for step size.
   I use White's values, although they seem very conservative, at least in terms of
   the sensitivity of the final answers |#

(define infty 5.)
(define h .02)

; We need a good guess at F'', which we find through interpolation
(define (F2->F1 F2)
  (let loop[(Y 0)
            (F 0)
            (F1 0)
            (F2 F2)]
    (if (< Y infty)
        (let-values([(F F1 F2)(integrate-F h F F1 F2)])
          (loop (+ Y h) F F1 F2))
        F1)))
(define small-F2 1.22)
(define small-F1 (F2->F1 small-F2))
(define big-F2 1.24)
(define big-F1 (F2->F1 big-F2))
(define F2 (+ small-F2 (/(*(- 1. small-F1)(- big-F2 small-F2))(- big-F1 small-F1))))
            
(module+ test
  (printf "Y\tF\tF1\tF2~n")
  (define h .1)
  (define (4digits x)
    (/ (round (* 1000. x)) 1000.))
  (let loop[(Y 0)
            (F 0)
            (F1 0)
            (F2 F2)]
    (printf"~a\t~a\t~a\t~a~n" (4digits Y)(4digits F)(4digits F1)(4digits F2))
    (cond[(< Y 10.)
          (let-values([(F F1 F2)(integrate-F h F F1 F2)])
            (loop (+ Y h) F F1 F2))])))

; refine F2
(set! small-F2 (- F2 .001))
(set! small-F1 (F2->F1 small-F2))
(set! big-F2 (+ F2 .001))
(set! big-F1 (F2->F1 big-F2))
(set! F2 (+ small-F2 (/(*(- 1. small-F1)(- big-F2 small-F2))(- big-F1 small-F1))))

; We accept the resulting solution as definitive
(define-values (Ys Fs F1s F2s)
  (let[(Y 0.)
       (F 0.)
       (F1 0.)
       (F2 F2)] 
    (for/lists
        (Ys Fs F1s F2s)
      [(count(in-range (/ infty h)))]
       ; to capture first point
      (let[(old-Y Y) 
           (old-F F)
           (old-F1 F1)
           (old-F2 F2)]
        (set!-values (F F1 F2)(integrate-F h F F1 F2))
        (set! Y (+ h Y))
        (values old-Y old-F old-F1 old-F2)))))

; Displacement thickness
(define D*0 (- (last Ys)(last Fs)))
(define (delta*0 B nu)
  (* D*0(sqrt (/ B nu))))

#| Momentum thickness
   Non-D version is integral wrt Y over interval [0, inf] of F'-(F')^2,
   if G'   = F'(1-F')
      G''  = F''(1-F')-F'F''=F''-2F'F''
      G''' = F'''-2(F'')^2-2F'F'''
   This sets us up well for my integral routine. |#
(define THETA0
  (- D*0
     (for/fold[(G 0.)]
              [(F Fs)
               (F1 F1s)
               (F2 F2s)]
       (let[(F3(F3 F F1 F2))]
         (integrate h G (* F1(- 1. F1))
                        (- F2(* 2. F1 F2))
                        (- F3 (* 2. (sqr F2))(* 2. F1 F3)))))))
       
(define (theta0 B nu)
  (* THETA0 (sqrt (/ B nu))))

; Starting shape factor
(define H0 (/ D*0 THETA0))

