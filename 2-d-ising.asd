(defsystem "2-d-ising"
  :description "2 维 Ising 模型的 Monte Carlo 模拟."
  :version "0.0.1"
  :author "凉凉 <l>"
  :licence "MIT"
  :depends-on ("gsll" "eazy-gnuplot" "lparallel")
  :components ((:file "mcmc-helper")
               (:file "metropolis" :depends-on ("mcmc-helper"))
               (:file "2-d-ising"  :depends-on ("mcmc-helper" "metropolis"))))
