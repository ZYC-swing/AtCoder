import java.io.OutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.Scanner;

/**
 * Built using CHelper plug-in
 * Actual solution is at the top
 *
 * @author ZYCSwing
 */
public class Main {
    public static void main(String[] args) {
        InputStream inputStream = System.in;
        OutputStream outputStream = System.out;
        Scanner in = new Scanner(inputStream);
        PrintWriter out = new PrintWriter(outputStream);
        TaskC solver = new TaskC();
        solver.solve(1, in, out);
        out.close();
    }

    static class TaskC {
        public void solve(int testNumber, Scanner in, PrintWriter out) {
            int a, b, c, x, y;
            a = in.nextInt();
            b = in.nextInt();
            c = in.nextInt();
            x = in.nextInt();
            y = in.nextInt();

            if (a + b <= 2 * c) {
                out.println(x * a + y * b);
            } else {
                int res;
                if (x < y) {
                    res = 2 * x * c + (y - x) * Math.min(b, 2 * c);
                } else {
                    res = 2 * y * c + (x - y) * Math.min(a, 2 * c);
                }
                out.println(res);
            }
        }

    }
}

