%------------------------------------------
%Algoritmo de clusters com agregação (HCA)
%------------------------------------------
hold off;
clear;
clc;
%Mensagem inicial(GUI) para introdução de elementos pretendidos
prompt = {'Introduza o número de elementos:'};
dlgtitle = 'Input';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
c=str2double(answer(1)); %Numero introduzido de agrupamentos

Cor=[]; %Cria um array que vai ser preenchido com valores random para corresponder a cor
for QQ=1:c
Cor=[Cor;rand(1,3)]; %Preenche a array com 3 valores random em linha
end

legenda= [1:c];

A=[];
for E=1:c %Criação de clusters de acordo com input
    a=[Cor(E,:)];
    Z = randn(10,2)+[randi(10),randi(10)];
    A=[A,Z];
    [n,nx]=size(A);
    figure(2);
    plot(Z(:,1),Z(:,2),'o','Color',a);
    text(Z(:,1),Z(:,2),string(legenda(E)));
    hold on;
end

myA=my_linkage(A,c,nx); %linkagem
figure(1);
title('Dendrograma HCA');
xlabel('Clusters') ;
ylabel('Distância') ;
dendrogram(myA); %dendograma final

%-------------------------------
%Metodo Top-Down
%-------------------------------
aux=c-1; %valor percorrer matriz G
Correspond=[]; % Matriz correspondencias de clusters aglomerados
Gnova=[];   %Matriz dist top-down

for i=c-1:-1:2 %linhas da matriz G da ultima ate primeira
   
    if (myA(i,2)>c)
        Correspond=[Correspond;myA(i-1,1),myA(i-1,2),myA(i,2)];
    end
    
    if (aux==c-1)
        Gnova=[Gnova;myA(aux,1),myA(aux,2),myA(aux,3)];
    end
    
    Gnova=[Gnova;myA(aux-1,1),myA(aux-1,2),myA(aux-1,3)];
    aux=aux-1;    
end
Gnova %Matriz Top-Down


%-------------------------------
%Metodo Bottom-Up
%-------------------------------
function G=my_linkage(A,c,nx)
%Matriz D vai possuir os valores de distancia entre agrupamentos
D=Inf(c); %Distancia de um cluster para ele proprio tem de ficar inutilizavel em contas
 for i=1:c
     for j=i+1:c

        cc_1=-1; %cluster comparacao coluna 1
        cc_2=0; %cluster comparacao coluna 2
        
        for l=1:c %ciclo for grande para percorrer quadrados todos
            p=1;
            cc_1=cc_1+2; 
            cc_2=cc_2+2;
            aux_1=1; %variavel auxiliar para percorrer a comparacao entre todas as colunas
            aux_2=2; %variavel auxiliar para percorrer a comparacao entre todas as colunas
            
            for p=1:c %ciclo for (pequeno) para saltar de linha 
             if cc_1==aux_1 && cc_2==aux_2
                 D(l,p)=Inf;
                 aux_1=aux_1+2;
                 aux_2=aux_2+2;
                 
             else
             D(l,p)=hausdorff(A(1:10,cc_1:cc_2),A(1:10,aux_1:aux_2)); %Calculo da distancia entre clusters  
             D(p,l)=D(l,p); %Distancia entre dois clusters é simétrica
             
             D; %Debug
            aux_1=aux_1+2;
            aux_2=aux_2+2;
             end
            end
        end
        
    end 
    end


G=[];
i=1;


%Minima dist entre clusters
while sum(sum(D<Inf))>1
    [m,I,J]=minMat(D);
    G=[G;min(I,J) max(I,J) m];
    minIJ=min(D([I,J],:));
    G
    D=[D minIJ'];
    D=[D;[minIJ Inf]];
    D(I,:)=Inf;
    D(J,:)=Inf;
    D(:,I)=Inf;
    D(:,J)=Inf;
    i=i+1;
end
end


function G=agrupamento_hier(X)

[n,nx]=size(X); % n. numero de pontos; nx - dimensão do espaço dos pontos

end

function [m,I,J]=minMat(D)
% Calcula o minimo elemento de D, indicando a sua posiçã linha (I) e coluna
% (J)
[mC,II]=min(D);
[m,J]=min(mC);
J=J(1);
I=II(J);
end

%-------------------------------
%Distancia de Hausdorff
%-------------------------------
%Calculo maximo entre distancias
function [dist] = hausdorff(A,B) 
dist = max(hdist(A,B), hdist(B,A));
end

%Cálculo mínimo entre distancias
function[dist] = hdist(A,B) 
MD=[];  
m = size(A, 1); 
n = size(B, 1); 
dim= size(A, 2); 
for k = 1:m 
    C = ones(n,1) * A(k,:); 
    Di = (C-B).*(C-B); 
    Di = sqrt(Di * ones(dim,1)); 
    dist(k) = min(Di);
    MD=[MD;dist(k)];
end
dist = max(MD(:,1));
end
