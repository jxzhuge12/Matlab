import java.util.Scanner;

public class bicycle {	
	public static void main(String[] args){
		Scanner in = new Scanner(System.in);
		int go = 1;
		do{
			System.out.println("输入满架数量：");
			int total;
			total = in.nextInt();//获取满架数量
			
			int i;
			int[] data = new int[100];
			
			System.out.println("输入车辆数据：");
			for(i = 0;i < 91;i++){
				data[i] = in.nextInt();
			}//获取91个时间点的数量
			
			int h = 0;
			for(i = 0;i < 91;i++){
				if(data[i] == 0){
					h = 1;
					int k = 0;
					int i2;
					for(i2 = i;i2 < 91;i2++)
						if(data[i2] != 0)break;//找出i后第一个不为0的数的位置，放入i2中
					for(int j = i;j < 91;j++){
						if(data[j] == total){
							k = 1;
							int lmax = 0;
							for(int t = i + 1;t < j;t++){
								int tt,ttt;
								for(tt = t;tt < j;tt++)
									if(data[tt] != data[t])break;
								for(ttt = t;ttt > i;ttt--)
									if(data[ttt] != data[t])break;
								if(data[t] > data[tt] && data[t] > data[ttt] && data[t] > lmax)
									lmax = data[t];
							}
							System.out.println("在 " + (6 + i / 6) + ":" + (i % 6) + "0 时，增加 " + Math.max(1,Math.min(Math.round((double)(total) / 4), (total - lmax - 1))) + " 辆车");
							for(int t = i;t < j;t++){
								data[t] += Math.max(1,Math.min(Math.round((double)(total) / 4), (total - lmax - 1)));
								if(data[t] > total)data[t] = total;
							}
							for(i = j;i > 0;i--){
								if(data[i] != total)break;
							}
							break;
						}
					}
					if(k == 0){
						int lmax = 0;
						for(int t = i + 1;t < 91;t++){
							int tt,ttt;
							for(tt = t;tt < 91;tt++)
								if(data[tt] != data[t])break;
							for(ttt = t;ttt > i;ttt--)
								if(data[ttt] != data[t])break;
							if(data[t] > data[tt] && data[t] > data[ttt] && data[t] > lmax)
								lmax = data[t];
						}
						if(data[90] > lmax)
							lmax = data[90];
						System.out.println("在 " + (6 + i / 6) + ":" + (i % 6) + "0 时，增加 " + Math.max(1,Math.min(Math.round((double)(total) / 4), (total - lmax - 1))) + " 辆车");
						for(int t = i;t < 91;t++){
							data[t] += Math.max(1,Math.min(Math.round((double)(total) / 4), (total - lmax - 1)));
							if(data[t] > total)
								data[t] = total;
						}
						i = i2 - 1;
					}//后来既没有0也没有最大值的情况
				}//遇到数据为0的情况
				
				else if(data[i] == total){
					h = 1;
					int k = 0;
					int i2;
					for(i2 = i;i2 < 91;i2++)
						if(data[i2] != total)break;//找出i后第一个不为total的数的位置，放入i2中
					for(int j = i;j < 91;j++){
						if(data[j] == 0){
							k = 1;
							int lmin = total;
							for(int t = i + 1;t < j;t++){
								int tt,ttt;
								for(tt = t;tt < j;tt++)
									if(data[tt] != data[t])break;
								for(ttt = t;ttt > i;ttt--)
									if(data[ttt] != data[t])break;
								if(data[t] < data[tt] && data[t] < data[ttt] && data[t] < lmin)
									lmin = data[t];
							}
							System.out.println("在 " + (6 + i / 6) + ":" + (i % 6) + "0 时，减少 " + Math.max(1,Math.min(Math.round((double)(total) / 4), (lmin - 1))) + " 辆车");
							for(int t = i;t < j;t++){
								data[t] -= Math.max(1,Math.min(Math.round((double)(total) / 4), (lmin - 1)));
								if(data[t] < 0)data[t] = 0;
							}
							for(i = j;i > 0;i--){
								if(data[i] != 0)break;
							}
							break;
						}//遇到数据为最大值后，后来出现0的情况
					}
					if(k == 0){
						int lmin = total;
						for(int t = i + 1;t < 91;t++){
							int tt,ttt;
							for(tt = t;tt < 91;tt++)
								if(data[tt] != data[t])break;
							for(ttt = t;ttt > i;ttt--)
								if(data[ttt] != data[t])break;
							if(data[t] < data[tt] && data[t] < data[ttt] && data[t] < lmin)
								lmin = data[t];
						}
						if(data[90] < lmin)
							lmin = data[90];
						System.out.println("在 " + (6 + i / 6) + ":" + (i % 6) + "0 时，减少 " + Math.max(1,Math.min(Math.round((double)(total) / 4), (lmin - 1))) + " 辆车");
						for(int t = i;t < 91;t++){
							data[t] -= Math.max(1,Math.min(Math.round((double)(total) / 4), (lmin - 1)));
							if(data[t] < 0)
								data[t] = 0;
						}
						i = i2 - 1;
					}//后来既没有0也没有最大值的情况
				}
			}
			if(h == 0){
				System.out.println("不需要对数据进行处理");
			}
		
			System.out.println("调整后的数据为:");
			for(i = 0;i < 90;i++){
				System.out.print(data[i] + ", ");
			}
			System.out.print(data[90]);
			System.out.println();
			
			System.out.println("1:重复程序;");
			System.out.println("0:退出程序.");
			go = in.nextInt();
		}while(go != 0);
	}
}
