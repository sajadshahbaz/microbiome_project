#######################################################################################################################
#Funkcja przelicza ciagla wartosc predykcji (np. prawdopodobienstwo) na decyzje binarna z zachowaniem proporcji klas
#v.1.0.0
#K. Mnich, 05.06.2019
#
#Parametry wejsciowe:
# x           wektor przewidywanego prawdopodobienstwa, ze zmienna decyzyjna bedzie 1 (np. randomforest$votes[,2])  
# y           wektor zmiennej decyzyjnej o wartosciach 0,1
#
#Wartosc funkcji - lista:
# threshold   wartosc poziomu odciecia
# decision    przewidywana wartosc zmiennej decyzyjnej
#
#
#Przyklad uzycia:
#
# model <- randomForest(x.training, y)
# bd <- balanced.decision(model$votes[,2], y)
# threshold <- bd$threshold
# predicted.decision.training <- bd$decision
#
# predicted.prob.validation <- predict(model, newdata=x.validation, type="prob")[,2]
# predicted.decision.validation <- predicted.prob.validation>=threshold
# cheated.decision.validation <- balanced.decision(predicted.prob.validation, y.validation)$decision
#
#UWAGA: 
# funkcja nie jest odporna na wystepowanie wielu jednakowych lub niemal jednakowych wartosci prawdopodobienstw
# np. jezeli wynik przewidywania bedzie albo bliski 0 albo bliski 1, to podzial moze wypasc w dziwnym miejscu
#######################################################################################################################

balanced.decision<-function(x, y) {
 
 if (any(y!=0 & y!=1)) stop('y should be 0 or 1') 
 if (length(x)!=length(y)) stop('x and y must have the same length')
 
 threshold<-x[order(x)[sum(y==0)]]
 
 decision<-x>=threshold
 
 return(list(threshold=threshold, decision=decision))
} 