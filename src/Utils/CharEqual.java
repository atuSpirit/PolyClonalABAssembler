package Utils;

public class CharEqual {
    public static boolean charEqual(char c1, char c2) {
        if (c1 == c2) {
            return true;
        }
        if (((c1 == 'I') && (c2 == 'L'))
                || ((c1 == 'L') && (c2 == 'I'))
                || ((c1 == 'N') && (c2 == 'D'))
                || ((c1 == 'D') && (c2 == 'N'))) {
            /*
                || ((c1 == 'Q') && (c2 == 'E'))
                || ((c1 == 'E') && (c2 == 'Q'))) {
                */

            return true;
        }
        return false;
    }
}
