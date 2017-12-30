#lang racket/base
;;; A collection of airfoil coordinates expressed as complex numbers
(provide (all-defined-out))
(module+ test
  (printf"airfoils.rkt running tests.~n")
  (require rackunit))
;;; I don't use regular expressions much, so breaking them out and testing them pays
(define empty-line-pattern  #px"^\\s*$")
(define 2-numbers-pattern   #px"(-?\\d*\\.\\d*)\\s+(-?\\d*\\.\\d*)$")
(define 1st-number-pattern  #px"^\\s*(-?\\d*\\.\\d*)")
(define lst-number-pattern  #px"(-?\\d*\\.\\d*)\\s*$")
(module+ test
  (check-true(regexp-match? empty-line-pattern "  "))
  (check-true(regexp-match? empty-line-pattern "\t \t"))
  (check-true(regexp-match? empty-line-pattern ""))
  (check-true(not(regexp-match? empty-line-pattern "hi there!")))
  (check-true(regexp-match? 2-numbers-pattern "-1. 2."))
  (check-true(regexp-match? 2-numbers-pattern ".2 -.1"))
  (check-true(regexp-match? 2-numbers-pattern "1.2 3.4"))
  (check-true(not(regexp-match? 2-numbers-pattern "1.2 Oo")))
  (check-true(regexp-match? 1st-number-pattern " .2 Ii "))
  (check-true(regexp-match? lst-number-pattern "Ee 1.")))


#| Translate *.dat files, such as those found at http://www.ae.illinois.edu/m-selig/ads/coord_database.html,
   to vectors of complex numbers for use by zfoil

   *.dat files ALWAYS have a line of text that usually has some textual description of the airfoil.
     "     "   ALWAYS consist only of numbers containing a decimal points after the first line.
   The second line ALWAYS consists of two numbers.
   IF the third line is blank, THEN the second line's numbers define the number of upper and lower surface points
                                    which will be listed from LE to TE, with upper and lower surfaces seperated 
                                    by a blank line.
                               ELSE the second and third lines are a part of a series of number pairs describing
                                    the surface CCW from and to the TE.
   The .dat files I've seen seem to favor FORTRAN conventions... there's always a decimal point somewhere. |#
(define(read-dat-file file-name)
  (define in (open-input-file file-name #:mode 'text))
  (define info-line(read-line in 'any))
  ; There are two varieties of .dat files, determined by nature of next two lines.
  (define-values(m n)(read-number-pair (read-line in 'any)))
  (define critical-line(read-line in 'any))
  (define zs
    (if(regexp-match? empty-line-pattern critical-line)
       ; m and n contains the number of upper and lower surface points.
       ; We expect to find upper and lower surface points, starting from LE, seperated by a blank line.
       (let[(upper(for/list[(i(in-range m))]
                    (let[(string(read-line in 'any))]
                      (let-values[((x y)(read-number-pair string))]
                        (make-rectangular x y)))))
            (blank(read-line in 'any)); It had better be blank!
            (lower(for/list[(i(in-range n))]
                    (let[(string(read-line in 'any))]
                      (let-values[((x y)(read-number-pair string))]
                        (make-rectangular x y)))))]
         (list->vector(append(reverse (cdr upper))lower)));cdr removes duplicate LE point
       ; If critical-line isn't empty, it had better hold another pair of numbers.
       ; Unknown numbers are the first coordinates of a list that continues until eof, defining
       ; points CCW from TE, just like zfoil!
       (list->vector(append
                     (list
                      (make-rectangular m n)
                      (let-values[((x y)(read-number-pair critical-line))]
                        (make-rectangular x y)))
                     (for/list[(l(in-lines in))]
                       (let-values[((x y)(read-number-pair l))]
                         (make-rectangular x y)))))))
  (close-input-port in)
  (values zs info-line))

;; read a pair of numbers
;; a FORTRAN artifact?
(define (read-number-pair string)
  (define strings(regexp-match 2-numbers-pattern string))
  (values(string->number(cadr strings))
         (string->number(caddr strings))))

(module+ test
  (let-values[((m n)(read-number-pair " 1.0 -2.0"))]
    (check-= m 1. 1e-99)
    (check-= n -2. 1e-99))
  #| These are not robust tests. They just demonstrate that the function runs without throwing
     an error. The output must be examined to confirm that things are as they should be.
     r1046 is of the first file type, that defines upper and lower surfaces seperately. e205 is
     of the second type, that defines points CCW from the trailing edge. |#
  (read-dat-file "Dat/r1046.dat")
  (read-dat-file "Dat/e205.dat"))

;; source: http://www.worldofkrauss.com/foils/793.dat
(define naca-4412 #(1.000000-0.000000i 0.960400+0.010627i 0.921600+0.020454i 0.883600+0.029531i
                    0.846400+0.037904i 0.810000+0.045612i 0.774400+0.052692i 0.739600+0.059175i
                    0.705600+0.065088i 0.672400+0.070456i 0.640000+0.075302i 0.608400+0.079646i
                    0.577600+0.083506i 0.547600+0.086901i 0.518400+0.089846i 0.490000+0.092358i
                    0.462400+0.094450i 0.435600+0.096139i 0.409600+0.097437i 0.384400+0.098326i
                    0.360000+0.098700i 0.336400+0.098576i 0.313600+0.097985i 0.291600+0.096957i
                    0.270400+0.095522i 0.250000+0.093709i 0.230400+0.091547i 0.211600+0.089066i
                    0.193600+0.086293i 0.176400+0.083257i 0.160000+0.079986i 0.144400+0.076506i
                    0.129600+0.072843i 0.115600+0.069024i 0.102400+0.065073i 0.090000+0.061014i
                    0.078400+0.056871i 0.067600+0.052665i 0.057600+0.048417i 0.048400+0.044148i
                    0.040000+0.039875i 0.032400+0.035616i 0.025600+0.031387i 0.019600+0.027202i
                    0.014400+0.023073i 0.010000+0.019012i 0.006400+0.015028i 0.003600+0.011130i
                    0.001600+0.007323i 0.000400+0.003612i 0.000000+0.000000i 
                    0.000400-0.003453i 0.001600-0.006685i 0.003600-0.009697i 0.006400-0.012489i
                    0.010000-0.015062i 0.014400-0.017417i 0.019600-0.019554i 0.025600-0.021475i
                    0.032400-0.023181i 0.040000-0.024675i 0.048400-0.025959i 0.057600-0.027036i
                    0.067600-0.027910i 0.078400-0.028584i 0.090000-0.029064i 0.102400-0.029356i
                    0.115600-0.029466i 0.129600-0.029401i 0.144400-0.029171i 0.160000-0.028786i
                    0.176400-0.028256i 0.193600-0.027594i 0.211600-0.026813i 0.230400-0.025929i
                    0.250000-0.024959i 0.270400-0.023920i 0.291600-0.022832i 0.313600-0.021718i
                    0.336400-0.020599i 0.360000-0.019500i 0.384400-0.018448i 0.409600-0.017457i
                    0.435600-0.016420i 0.462400-0.015316i 0.490000-0.014158i 0.518400-0.012962i
                    0.547600-0.011742i 0.577600-0.010516i 0.608400-0.009297i 0.640000-0.008102i
                    0.672400-0.006945i 0.705600-0.005842i 0.739600-0.004804i 0.774400-0.003843i
                    0.810000-0.002968i 0.846400-0.002186i 0.883600-0.001502i 0.921600-0.000913i
                    0.960400-0.000416i 1.000000+0.000000i))
;; source: http://wind.nrel.gov/public/library/3387.pdf
(define naca-66_2-415 #(1.00000+0.00000i 0.95053+0.01196i 0.90104+0.02519i 0.85139+0.03872i
                        0.80159+0.05187i 0.75162+0.06419i 0.70150+0.07518i 0.65126+0.08431i
                        0.60090+0.09100i 0.55046+0.09473i 0.50000+0.09656i 0.44952+0.09685i
                        0.39904+0.09571i 0.34857+0.09309i 0.29812+0.08897i 0.24771+0.08329i
                        0.19736+0.07581i 0.14709+0.06624i 0.09696+0.05381i 0.07199+0.04617i
                        0.04711+0.03718i 0.02241+0.02592i 0.01019+0.01873i 0.00544+0.01467i
                        0.00314+0.01206i 0.00000+0.00000i 0.00686-0.01006i 0.00956-0.01187i
                        0.01481-0.01445i 0.02759-0.01848i 0.05289-0.02454i 0.07801-0.02921i
                        0.10304-0.03313i 0.15291-0.03932i 0.20264-0.04397i 0.25229-0.04749i
                        0.30188-0.05009i 0.35143-0.05189i 0.40096-0.05287i 0.45048-0.05305i
                        0.50000-0.05244i 0.54954-0.05093i 0.59910-0.04816i 0.64874-0.04311i
                        0.69850-0.03630i 0.74838-0.02839i 0.79841-0.02003i 0.84861-0.01180i
                        0.89896-0.00451i 0.94947+0.00068i 1.00000+0.00000i))
;; Generate a Joukowski airfoil
(require racket/math
         "Private/Inviscid/joukowsky.rkt")
(define (joukowsky-foil j0 c0 wte n)
  (define-values (zte a b)(make-dependent-parameters c0 j0 wte))
  (define dPhi (/(* 2 pi)n))
  (define C (make-circle c0 a))
  (define J (JofZ j0 b))
  (build-vector (add1 n) (λ(i) (J(C(* dPhi i))))))
(module+ test
  (require "Private/Utility/plots.rkt")
  (plot-foil (joukowsky-foil 0.5 0.52+0.01i 1. 21)))
(module+ main
  (require(submod ".." test)))