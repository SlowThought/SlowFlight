#lang racket/base
(require racket/math rackunit)
(provide powell test-powell)
;; Numerical optimization routines for use in zfoil
(define epsilon 1e-14)
;; Test function is meant to be somewhat like airfoil-fit-error, with "valley" skewed wrt x, y axes
;; By inspection, minimum at (0 0).

(define (test-function x y)
  (+ (sqr x)
     (sqr y)
     (* x y)))
(require plot)
(define (test-test-function)
  (display(plot3d (surface3d test-function)
                #:x-min -1.0
                #:x-max  1.0
                #:y-min -1.0
                #:y-max  1.0))
  (printf "~nYou should see a pretty straightforward minimization problem~n"))

;; Find minimum of quadratic function f(t) fitted to t = -1, 0, 1 using determinants
(define (find-minimum fl fm fr)
  (/(- fl fr)
    (- (+ fr fl)(* 2 fm))
    2))
(define (test-find-minimum)
  (check-= (find-minimum 1 0 1) 0 epsilon)
  (check-= (find-minimum 1 4 9) -2 epsilon)
  (check-= (find-minimum 4 16 36) -2 epsilon))

;; Find minimum of real function of a complex variable
;; Use of parallelism improves speed, but f better not have side effects!
(require racket/future)
(define (powell f z0 h epsilon)
  ; First we move in h direction
  (let*[(f1(future(λ()(f(- z0 h)))))
        (f2(future(λ()(f z0))))
        (z1(+ z0(* h(find-minimum(touch f1)
                                 (touch f2)
                                 (f(+ z0 h))))))]
    ; Then we move normal to it
    (let*[(i(* 0+i h))
          (f3(future(λ()(f(- z1 i)))))
          (f4(future(λ()(f z1))))
          (z2(+ z1(* i(find-minimum(touch f3)
                                   (touch f4)
                                   (f(+ z1 i))))))
          ; Then we find a new step that sort of is near the steepest gradient
          (new-h(- z2 z0))]
      (if (< (magnitude new-h) epsilon); we're done
          z2; and we return the point corresponding to the minimum
          (powell f z2 new-h epsilon))))); else we search some more in a hopefully better direction
(define (test-powell-function)
  (check-= (magnitude(powell (λ(z)(test-function(real-part z)(imag-part z)))-10.0+100.0i 1.0 epsilon))
           0.0 epsilon))
;; Test code
(define (test-powell)
  (test-test-function)
  (test-powell-function))

;; Parallel test code
(define (expensive-function z)
  (for[(i 100)]
    (test-function(real-part z)(imag-part z)))
  (test-function(real-part z)(imag-part z)))
(define (time-powell)
  (time
   (for[(i 10)]
     (powell expensive-function -10.0+100.0i 1.0 epsilon))))

  
