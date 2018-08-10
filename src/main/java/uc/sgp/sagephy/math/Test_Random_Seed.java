package uc.sgp.sagephy.math;

import java.util.Random;
import java.util.concurrent.locks.ReentrantLock;
import org.uncommons.maths.binary.BinaryUtils;
import org.uncommons.maths.random.RepeatableRNG;
import org.uncommons.maths.random.SeedException;
import org.uncommons.maths.random.SeedGenerator;

public class Test_Random_Seed {

	public static void main(String[] args) {
		
		for(int i = 0; i < 100; i++) {
			byte[] seed = DefaultSeedGenerator2.getInstance().generateSeed(16);
			System.out.println(i + ": " + seed);
		}

	}

}
