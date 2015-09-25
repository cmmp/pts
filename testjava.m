%load_javaplex;

%ari = MyUtils.computeAdjustedRandIndex([1 0 1], [1 1 2]);

%ari



X = csvread('C:/Users/Cássio/workspace/darbseg/resources/Y2.csv');


javaaddpath('C:\Users\Cássio\workspace\gridseg\target\gridseg-0.0.1-SNAPSHOT-jar-with-dependencies.jar');
import gridseg.*;

lbls = gridseg.GridSeg.gridseg(X, 35);


