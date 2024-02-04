function dy=Genericlorenz(t,y)
dy=zeros(3,1);
dy(1)=10*(-y(1)+y(2));
dy(2)=28*y(1)-y(2)-1*y(1)*y(3);
dy(3)=1*(y(1)*y(2)-(8/3)*y(3));
end
