function W_matrix = calculate_W_matrix(i,Db,Dr)
% функция для расчета матрицы W
    M1 = [Db(1,i) 0 Db(3,i) -Db(2,i);
          0 Db(1,i) Db(2,i) Db(3,i);
          Db(3,i) Db(2,i) -Db(1,i) 0;
          -Db(2,i) Db(3,i) 0 -Db(1,i)];
    M2 = [Db(2,i) -Db(1,i) 0 Db(1,i);
          -Db(3,i) -Db(2,i) Db(1,i) 0;
          0 Db(1,i) Db(2,i) Db(3,i);
          Db(1,i) 0 Db(3,i) -Db(2,i)];
    M3 = [Db(3,i) Db(2,i) -Db(1,i) 0;
          Db(2,i) -Db(3,i) 0 Db(1,i);
          -Db(1,i) 0 -Db(3,i) Db(2,i);
          0 Db(1,i) Db(2,i) Db(3,i)];
    W_matrix = Dr(1,i)*M1 + Dr(2,i)*M2 + Dr(3,i)*M3;
end