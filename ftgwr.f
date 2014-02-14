program ftgwr
	implicit none
	integer, parameter :: prec = 8
	integer, parameter :: prec_int = 8
	complex (kind=prec) :: t, p, x, aux, an, bn
	real (kind=prec) :: u, k, z, q, h
	integer (kind=prec_int) :: i, M
	x = (5.5, 0)
	t = (2.0, 0)

	!roda para vários M diferentes para ver a precisão, em geral com M=32 está OK
   	M = 1
   	!pause
	do i=M, 50
		write(*,*) 'M=',i,' f(x,t)=',gwr(x, t, i)
		!pause
	end do
	
	contains
	
	!Função que se deseja inverter duas vezes
	complex (kind=prec) function func(p, s)
		complex (kind=prec), intent(in) :: p, s
		!func = -Log(p*s)/(p*s)
		func = 1./(s*s * p*p)
		!func = 1./(s*p)
	end function func

	complex (kind=prec) function Sm(theta, r)
		implicit none
		real (kind=prec), intent(in) :: theta
		complex (kind=prec), intent(in) :: r
		complex (kind=prec) :: temp
		Sm = r*theta*complex((1.0/tan(theta)),1.0)
	end function Sm
	
	real (kind=prec) function Si(theta)
		implicit none
		real (kind=prec), intent(in) :: theta
		Si = theta + (theta * (1./tan(theta)) - 1) * (1./tan(theta))
	end function Si
	
	complex (kind=prec) function FT(t, s, M)
		implicit none
		complex (kind=prec), intent(in) :: t, s
		complex (kind=prec) :: r
		real (kind=prec) :: su, thetak, pi
		integer (kind=prec_int) :: k
		integer (kind=prec_int), intent(in) :: M
		
		pi = dacos(DBLE(-1.))
		r = (2*M) / (301*t)
		
		su = 0.
		do k=1,M-1
			thetak = k*pi/M
			su = su + real(exp(t * Sm(thetak, r)) * func(Sm(thetak, r), s) * complex(1.0, Si(thetak)))
		end do
		
		FT = (r/M) * (func(r, s)/2. * exp(r*t) + su)
	end function FT

	!O GWR2 só deve ser utilizado pelo GWR, para inverter a primeira variável
	!Mostrou precisão ruim, deve ser pelo 4*M
	complex (kind=prec) function gwr2(s, p, M)
		implicit none
		complex (kind=prec):: expr, sum, s1, best
		complex (kind=prec), intent(in) :: s, p
		complex (kind=prec), dimension(2*M) :: Fi
		complex (kind=prec), dimension(M) :: G, Gp, Gm
		integer (kind=prec_int) :: k, n, i, M1
		integer (kind=prec_int), intent(in) :: M
		real (kind=prec) :: temp
		
		s1 = 0.6931471805599453094172321214581765680755001343602552541206800094934/s
		
		do i=1,2*M
			Fi(i) = func(i*s1, p)	
			!print*, Fi(i)
		end do
		
		M1 = M
		
		do n=1,M
			sum = 0
			do i = 0, n
				sum = sum + (binomial(n, i) * ((-1.) ** i) * Fi(n+i))
			end do
			
			G(n) = s1 * sum * (fat_int(2*n) / (fat_int(n)*fat_int(n-1)))
			!print*, (fat_int(2*n) / (fat_int(n)*fat_int(n-1)))
		end do
		
		do n = 0, M1
			Gm(n+1) = 0
		end do
		best = G(M1)
		
		do k=0,M1-2
			do n = M1-2-k, 0, -1
				expr = G(n+2) - G(n+1)
				if (expr == (0., 0.)) then
					gwr2 = best
					return
				endif
				expr = Gm(n+2) + (k + 1)/expr
				Gp(n+1) = expr
				
				if (mod(k, 2) == 1) then
					if (n == M1 - 2 - k) then
						best = expr
					end if
				end if
			end do
			
			do n=0,M1-k
				Gm(n+1) = G(n+1)
				G(n+1) = Gp(n+1)
			end do
		end do

		gwr2 = best
	end function gwr2
	
	!Inverte a segunda variável
	complex (kind=prec) function gwr(s, p, M)
		implicit none
		complex (kind=prec):: expr, sum, p1, best
		complex (kind=prec), intent(in) :: p, s
		complex (kind=prec), dimension(2*M) :: Fi
		complex (kind=prec), dimension(M) :: G, Gp, Gm
		integer (kind=prec_int) :: k, n, i, M1
		integer (kind=prec_int), intent(in) :: M
		real (kind=prec) :: temp
		
		p1 = 0.6931471805599453094172321214581765680755001343602552541206800094934/p
		
		!chama FT ou GWR2 para a primeira variável
		do i=1,2*M
			!Fi(i) = gwr2(s, i*p1, 4*M)		!Chama GWR
			Fi(i) = FT(s, i*p1, 3*M)		!Chama FT
			!print*, real(Fi(i))
			!pause
		end do
		
		M1 = M;
		
		do n=1,M
			sum = 0
			do i = 0, n
				sum = sum + (binomial(n, i) * ((-1) ** i) * Fi(n+i))
			end do
			
			G(n) = p1 * sum * (fat_int(2*n) / (fat_int(n)*fat_int(n-1)))
			!write(*,*) 'g',real(G(n))
		end do
		
		do n = 0, M1
			Gm(n+1) = 0
		end do
		best = G(M1)
		
		do k=0,M1-2
			do n = M1-2-k, 0, -1
				expr = G(n+2) - G(n+1)
				if (expr == (0., 0.)) then
					gwr = best
					return
				endif
				expr = Gm(n+2) + (k + 1)/expr
				Gp(n+1) = expr
				
				if (mod(k, 2) == 1) then
					if (n == M1 - 2 - k) then
						best = expr
					end if
				end if
				
			end do
			
			do n=0,M1-k
				Gm(n+1) = G(n+1)
				G(n+1) = Gp(n+1)
			end do
		end do
	
		gwr = best
	end function gwr

	!Para fatorial de números reais
	real (kind=prec) function fat(x)
		implicit none
		real (kind=prec), intent(in) :: x
		fat = dgamma(DBLE(x+1.0))
	end function fat

	!Para fatorial de números inteiros, aumenta a precisão
	real (kind=prec) function fat_int(x)
		implicit none
		integer (kind=prec_int), intent(in) :: x
		integer (kind=prec_int) :: y
		real (kind=prec) :: cont
		y = x
		cont = 1
		do while (y > 1)
			cont = cont*y
			y = y - 1
		end do
		fat_int = cont
	end function fat_int

	real (kind=prec) function binomial(n, k)
		implicit none
		integer (kind=prec_int), intent(in) :: n,k
		binomial = fat_int(n)/(fat_int(k)*fat_int(n-k))
	end function binomial
	
end program
